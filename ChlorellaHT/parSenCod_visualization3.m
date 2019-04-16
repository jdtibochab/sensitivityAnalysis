%% Including lipid accumulation

clear
% load('parSenCod_PA_50BOF.mat')
% load('paralleltestPA.mat')
% load('paralleltestPA_WS.mat')
% load('parSenCod_PA.mat')
% load('paralleltestPA_CZ.mat')
% load('parSenCod_HT.mat')
% load('parSenCod_Pt.mat')
% load('parSenCod_PA_WS2.mat')
% load('paralleltestPt_WS.mat')
% load('paralleltestPA_complete.mat')
load('paralleltestPt.mat')

% model=model_an;
Met2=Met;
model=model_an;

% Lipids
TAG_mets = strmatch('tag',Met);
DAG_mets = strmatch('dgdg',Met);
MAG_mets = strmatch('mgdg',Met);
Pe_mets = strmatch('pe',Met);
Pg_mets = strmatch('pg',Met);
As_mets = strmatch('asqd',Met);
Dgts_mets = strmatch('dgts',Met);
Dgdg_mets = strmatch('dgdg',Met);
    Lip_mets = sort([TAG_mets;DAG_mets;MAG_mets;Pe_mets;Pg_mets;As_mets;Dgts_mets]);

for i=1:length(TAG_mets)
    Met2(TAG_mets(i))={['tag',num2str(i)]};
end
for i=1:length(DAG_mets)
    Met2(DAG_mets(i))={['dag',num2str(i)]};
end
for i=1:length(MAG_mets)
    Met2(MAG_mets(i))={['mag',num2str(i)]};
end
for i=1:length(Dgts_mets)
    Met2(Dgts_mets(i))={['dgts',num2str(i)]};
end
for i=1:length(Dgdg_mets)
    Met2(Dgdg_mets(i))={['dgdg',num2str(i)]};
end
for i=1:length(As_mets)
    Met2(As_mets(i))={['asqd',num2str(i)]};
end
for i=1:length(Pg_mets)
    Met2(Pg_mets(i))={['pg',num2str(i)]};
end
for i=1:length(Pe_mets)
    Met2(Pe_mets(i))={['pe',num2str(i)]};
end
    
% AA
AA_mets1 = strfind(Met,'__L');
AA_mets2 = strfind(Met,'gly_c');

j = 1;
for i = 1:length(AA_mets1)
   if AA_mets2{i} == 1
       AA_mets(j,1) = i;
       j = j+1;
   end
   if AA_mets1{i} == 4
       AA_mets(j,1) = i;
       j = j+1;
   end
end
for i = 1:length(AA_mets)
   Met2(AA_mets(i)) = strtok(Met(AA_mets(i)),'_');    
end
AA_mets(strmatch('gal__L_c',Met(AA_mets))) = [];

% CB
CB_mets1 = strfind(model.metNames(Met_ID),'ose');
notCB_mets1 = strfind(model.metNames(Met_ID),'serine');
j = 1;
for i = 1:length(CB_mets1)
   if CB_mets1{i} ~= 0
       if notCB_mets1{i} ~= 0
       else
           CB_mets(j,1) = i;
           j = j+1;
       end
   end
end
CB_mets(end) = []; % starch

% Nuc
Nuc_mets1 = strfind(Met,'tp');
j = 1;
for i = 1:length(Nuc_mets1)
   if Nuc_mets1{i} ~= 0
       Nuc_mets(j,1) = i;
       j = j+1;
   end
end
Nuc_mets(1) = []; % ATP

%% By groups
group = 1:length(Met);
% group = AA_mets; groupname = 'AA';
% group = TAG_mets; groupname = 'TAG';
% group = DAG_mets; groupname = 'DAG';
% group = MAG_mets; groupname = 'MAG';
% group = Pe_mets; groupname = 'Pe';
% group = Pg_mets; groupname = 'Pg';
% group = As_mets; groupname = 'As';
% group = Dgts_mets; groupname = 'Dgts';
% group = Dgdg_mets; groupname = 'Dgdg';
% group = CB_mets; groupname = 'CB';


for i = 1:length(group)
    k = group(i);
    co2_uptake{k} = abs(co2_uptake{k});
    parMetCoeff{k} = abs(parMetCoeff{k});
    for j=1:it1 % For all BOFs
        growth_mean(i,j) = (max(growth{k}(:,j))+min(growth{k}(:,j)))/2; % Mean value of all points of 1 BOF
        growth_change(i,j)= max(growth{k}(:,j))-min(growth{k}(:,j)); % Change within 1 BOF
        parCoeff_change(i,j) = max(parMetCoeff{k}(:,j))-min(parMetCoeff{k}(:,j)); % Change within 1 BOF
        co2_change(i,j) = max(co2_uptake{k}(:,j))-min(co2_uptake{k}(:,j)); % Change within 1 BOF
        [poly] = polyfit(parMetCoeff{k}(:,j),co2_uptake{k}(:,j),1); 
    end
    
    min_change = min(growth_change(i,:));
    max_change = max(growth_change(i,:));
    
    N_temp = find(growth_change(i,:)==max_change);
    N(i) = N_temp(1);
    
    N2_temp = find(growth_change(i,:)==min_change);
    N2(i) = N2_temp(1);
    
    growth_sens_max(i) = growth_change(i,N(i))/parCoeff_change(i,N(i));
%     plot(N(i),growth_sens_max(i),'^','Color',colormat(i,:),'MarkerFaceColor',colormat(i,:))
    
    
    growth_sens_min(i) = growth_change(i,N2(i))/parCoeff_change(i,N2(i));
%     plot(N2(i),growth_sens_min(i),'s','Color',colormat(i,:),'MarkerFaceColor',colormat(i,:))
end

%%
growth_sens = mean([growth_sens_min;growth_sens_max]);

AAsens = mean(growth_sens(AA_mets));
    AAsens_dev = std(growth_sens(AA_mets));
TAGsens = mean(growth_sens(TAG_mets));
    TAGsens_dev = std(growth_sens(TAG_mets));
DAGsens = mean(growth_sens(DAG_mets));
    DAGsens_dev = std(growth_sens(DAG_mets));
MAGsens = mean(growth_sens(MAG_mets));
    MAGsens_dev = std(growth_sens(MAG_mets));
Pesens = mean(growth_sens(Pe_mets));
    Pesens_dev = std(growth_sens(Pe_mets));
Pgsens = mean(growth_sens(Pg_mets));
    Pgsens_dev = std(growth_sens(Pg_mets));
Assens = mean(growth_sens(As_mets));
    Assens_dev = std(growth_sens(As_mets));
Dgtssens = mean(growth_sens(Dgts_mets));
    Dgtssens_dev = std(growth_sens(Dgts_mets));
Dgdgsens = mean(growth_sens(Dgdg_mets));
    Dgdgsens_dev = std(growth_sens(Dgdg_mets));
CBsens = mean(growth_sens(CB_mets));
    CBsens_dev = std(growth_sens(CB_mets));
Nucsens =  mean(growth_sens(Nuc_mets));
    Nucsens_dev = std(growth_sens(Nuc_mets));
ATPsens = growth_sens(strmatch('atp_c',Met));
    ATPsens_dev = 0;
    
Sens = [AAsens,TAGsens,DAGsens,MAGsens,Pesens,Pgsens,Assens,Dgtssens,Dgdgsens,CBsens,Nucsens];
Sens_dev = [AAsens_dev,TAGsens_dev,DAGsens_dev,MAGsens_dev,Pesens_dev,Pgsens_dev,Assens_dev,Dgtssens_dev,Dgdgsens_dev,CBsens_dev,Nucsens_dev];

%%
growth_cod = abs(growth_sens_min-growth_sens_max);

AAcod = mean(growth_cod(AA_mets));
    AAcod_dev = std(growth_cod(AA_mets));
TAGcod = mean(growth_cod(TAG_mets));
    TAGcod_dev = std(growth_cod(TAG_mets));
DAGcod = mean(growth_cod(DAG_mets));
    DAGcod_dev = std(growth_cod(DAG_mets));
MAGcod = mean(growth_cod(MAG_mets));
    MAGcod_dev = std(growth_cod(MAG_mets));
Pecod = mean(growth_cod(Pe_mets));
    Pecod_dev = std(growth_cod(Pe_mets));
Pgcod = mean(growth_cod(Pg_mets));
    Pgcod_dev = std(growth_cod(Pg_mets));
Ascod = mean(growth_cod(As_mets));
    Ascod_dev = std(growth_cod(As_mets));
Dgtscod = mean(growth_cod(Dgts_mets));
    Dgtscod_dev = std(growth_cod(Dgts_mets));
Dgdgcod = mean(growth_cod(Dgdg_mets));
    Dgdgcod_dev = std(growth_cod(Dgdg_mets));
CBcod = mean(growth_cod(CB_mets));
    CBcod_dev = std(growth_cod(CB_mets));
Nuccod = mean(growth_cod(Nuc_mets));
    Nuccod_dev = std(growth_cod(Nuc_mets));
ATPcod = growth_cod(strmatch('atp_c',Met));
    ATPcod_dev = 0;
    
Cod = [AAcod,TAGcod,DAGcod,MAGcod,Pecod,Pgcod,Ascod,Dgtscod,Dgdgcod,CBcod,Nuccod];
Cod_dev = [AAcod_dev,TAGcod_dev,DAGcod_dev,MAGcod_dev,Pecod_dev,Pgcod_dev,Ascod_dev,Dgtscod_dev,Dgdgcod_dev,CBcod_dev,Nuccod_dev];

%%

met_an = {'AA';'TAG';'DAG';'MAG';'PE';'PG';'AS';'DGTS';'DGDG';'CB';'N'};

%% Correlation for A_nutrient

groupCost = biosyntheticCost(model,Met(group),0.025,'EX_co2(e)'); % PA
% groupCost = biosyntheticCost(model,Met(group),0.025,'EX_glc__D_e'); % Pt
% groupCost = biosyntheticCost(model,Met(group),0.025,'EX_glc-A(e)'); % HT

AAcost = mean(groupCost(AA_mets));
    AAcost_dev = std(groupCost(AA_mets));
TAGcost = mean(groupCost(TAG_mets));
    TAGcost_dev = std(groupCost(TAG_mets));
DAGcost = mean(groupCost(DAG_mets));
    DAGcost_dev = std(groupCost(DAG_mets));
MAGcost = mean(groupCost(MAG_mets));
    MAGcost_dev = std(groupCost(MAG_mets));
Pecost = mean(groupCost(Pe_mets));
    Pecost_dev = std(groupCost(Pe_mets));
Pgcost = mean(groupCost(Pg_mets));
    Pgcost_dev = std(groupCost(Pg_mets));
Ascost = mean(groupCost(As_mets));
    Ascost_dev = std(groupCost(As_mets));
Dgtscost = mean(groupCost(Dgts_mets));
    Dgtscost_dev = std(groupCost(Dgts_mets));
Dgdgcost = mean(groupCost(Dgdg_mets));
    Dgdgcost_dev = std(groupCost(Dgdg_mets));
CBcost = mean(groupCost(CB_mets));
    CBcost_dev = std(groupCost(CB_mets));
Nuccost = mean(groupCost(Nuc_mets));
    Nuccost_dev = std(groupCost(Nuc_mets));
ATPcost = groupCost(strmatch('atp_c',Met));
    ATPcost_dev = 0;
    
Cost = [AAcost,TAGcost,DAGcost,MAGcost,Pecost,Pgcost,Ascost,Dgtscost,Dgdgcost,CBcost,Nuccost];
Cost_dev = [AAcost_dev,TAGcost_dev,DAGcost_dev,MAGcost_dev,Pecost_dev,Pgcost_dev,Ascost_dev,Dgtscost_dev,Dgdgcost_dev,CBcost_dev,Nuccost_dev];

figure;
subplot(1,2,1)
hold on
errorbar(Cost,Sens,Sens_dev,Sens_dev,Cost_dev,Cost_dev,'*')
title('Sensitivity')
text(Cost+0.5,Sens,met_an,'FontSize',8)
xlabel('A_{CO_2}')
ylabel('\Delta\mu/\Deltap')
xlim([0 70])
hold off

subplot(1,2,2)
hold on
errorbar(Cost,Cod,Cod_dev,Cod_dev,Cost_dev,Cost_dev,'*')
title('Parameter codependence')
text(Cost+0.5,Cod,met_an,'FontSize',8)
xlabel('A_{CO_2}')
ylabel('\Delta(\Delta\mu/\Deltap)')
x0=10;
y0=10;
width=1000;
height=800;
xlim([0 70])
hold off
set(gcf,'units','points','position',[x0,y0,width,height])

saveas(gcf,'summaryA','svg')

%% Correlation with conectivity
groupConectivity = calculateConectivity(model,Met(group));

AA_mets(strmatch('glu__L_c',Met(AA_mets))) = [];

met_an = {'AA';'glu__L_c';'TAG';'DAG';'MAG';'PE';'PG';'AS';'DGTS';'DGDG';'CB';'N'};

AAsens = mean(growth_sens(AA_mets));
    AAsens_dev = std(growth_sens(AA_mets));
Glusens = growth_sens(strmatch('glu__L_c',Met));
    Glusens_dev = 0;
AAcod = mean(growth_cod(AA_mets));
    AAcod_dev = std(growth_cod(AA_mets));
Glucod = growth_cod(strmatch('glu__L_c',Met));
    Glucod_dev = 0;

Sens = [AAsens,Glusens,TAGsens,DAGsens,MAGsens,Pesens,Pgsens,Assens,Dgtssens,Dgdgsens,CBsens,Nucsens];
    Sens_dev = [AAsens_dev,Glusens_dev,TAGsens_dev,DAGsens_dev,MAGsens_dev,Pesens_dev,Pgsens_dev,Assens_dev,Dgtssens_dev,Dgdgsens_dev,CBsens_dev,Nucsens_dev];
Cod = [AAcod,Glucod,TAGcod,DAGcod,MAGcod,Pecod,Pgcod,Ascod,Dgtscod,Dgdgcod,CBcod,Nuccod];
    Cod_dev = [AAcod_dev,Glucod_dev,TAGcod_dev,DAGcod_dev,MAGcod_dev,Pecod_dev,Pgcod_dev,Ascod_dev,Dgtscod_dev,Dgdgcod_dev,CBcod_dev,Nuccod_dev];


AAcon = mean(groupConectivity(AA_mets));
    AAcon_dev = std(groupConectivity(AA_mets));
Glucon = groupConectivity(strmatch('glu__L_c',Met));
    Glucon_dev = 0;
TAGcon = mean(groupConectivity(TAG_mets));
    TAGcon_dev = std(groupConectivity(TAG_mets));
DAGcon = mean(groupConectivity(DAG_mets));
    DAGcon_dev = std(groupConectivity(DAG_mets));
MAGcon = mean(groupConectivity(MAG_mets));
    MAGcon_dev = std(groupConectivity(MAG_mets));
Pecon = mean(groupConectivity(Pe_mets));
    Pecon_dev = std(groupConectivity(Pe_mets));
Pgcon = mean(groupConectivity(Pg_mets));
    Pgcon_dev = std(groupConectivity(Pg_mets));
Ascon = mean(groupConectivity(As_mets));
    Ascon_dev = std(groupConectivity(As_mets));
Dgtscon = mean(groupConectivity(Dgts_mets));
    Dgtscon_dev = std(groupConectivity(Dgts_mets));
Dgdgcon = mean(groupConectivity(Dgdg_mets));
    Dgdgcon_dev = std(groupConectivity(Dgdg_mets));
CBcon = mean(groupConectivity(CB_mets));
    CBcon_dev = std(groupConectivity(CB_mets));
Nuccon = mean(groupConectivity(Nuc_mets));
    Nuccon_dev = std(groupConectivity(Nuc_mets));

Con = [AAcon;Glucon;TAGcon;DAGcon;MAGcon;Pecon;Pgcon;Ascon;Dgtscon;Dgdgcon;CBcon;Nuccon]';
    Con_dev = [AAcon_dev,Glucon_dev,TAGcon_dev,DAGcon_dev,MAGcon_dev,Pecon_dev,Pgcon_dev,Ascon_dev,Dgtscon_dev,Dgdgcon_dev,CBcon_dev,Nuccon_dev];

figure;
subplot(1,2,1)
errorbar(Con,Sens,Sens_dev,Sens_dev,Con_dev,Con_dev,'*')
xlim([0,20])
title('Sensitivity')
text(Con+0.15,Sens+1.5e-4,met_an,'FontSize',8)
xlabel('Conectivity')
ylabel('\Delta\mu/\Deltap')


subplot(1,2,2)
errorbar(Con,Cod,Cod_dev,Cod_dev,Con_dev,Con_dev,'*')
xlim([0,20])
title('Parameter codependence')
text(Con+0.15,Cod+1.5e-4,met_an,'FontSize',8)
xlabel('Conectivity')
ylabel('\Delta(\Delta\mu/\Deltap)')
x0=10;
y0=10;
width=1000;
height=800;
set(gcf,'units','points','position',[x0,y0,width,height])

saveas(gcf,'summaryConec2','svg')