clear
% Pt
tic
fprintf('\n\n Calculating for Yl TEST1\n\n')

load('iYL_ready.mat')
modelYl = modelD;
modelYl = changeObjective(modelYl,'xBIOMASS');
% out = optimizeCbModel(modelD)

[Stoich,Met]=xlsread('bofYl_alt.xls');
Stoich(1,:)=[];
Met(1)=[];

% [model_an,prodCapMet,notProdCapMet,Stoech2] = prodCap(modelYl,Met,Stoich);
[model_an,prodCapMet,notProdCapMet,Stoich2,~] = metTest(modelYl,Met,Stoich);

Stoich = -[Stoich2;ones(1,length(Stoich(1,:)))*86.7881;ones(1,length(Stoich(1,:)))*86.7881];
prodCapMet = [prodCapMet;'atp_c';'h2o_c'];

Met = prodCapMet;
%%
t0 = [0.5 1.625 3.62];
% t = [0.5 1 1.5 2 3.5 5];
t = [min(t0):(max(t0)-min(t0))/5:max(t0)];
[intStoich,rsq] = interpolateStoich(Stoich,t0,t);
Stoich = intStoich;

%% Biomass composition change
mets = Met;
model = model_an;
for i = 1:length(Stoich(1,:))
    molarCoeffs = abs(Stoich(:,i));
    molarCoeffs(molarCoeffs>1) = 0;
    metIDs = findMetIDs(model,mets);
    [~,molarMass] = calculateFormula(model,metIDs,1);
    massCoeffs = molarCoeffs.*molarMass;
    totalMass(i) = sum(massCoeffs);
    massFrac{i} = massCoeffs/totalMass(i);
end
[metTypes,mets2] = findMetType(model,mets);

metGroups.CB = {'CB'};
metGroups.rest = {'rest'};
metGroups.Nuc = {'Nuc'};
metGroups.AA = {'AA'};
metGroups.MAG = {'MAG'};
metGroups.TAG = {'TAG'};
metGroups.PL = {'Pe','Pg','Pail','Pchol'};
metGroups.OL = {'As','Dgts','Dgdg','Sqdg','Mgdg'};
metFieldNames = fieldnames(metGroups);
for i = 1:length(Stoich(1,:))
    for m = 1:length(metFieldNames)
        fieldNames = metGroups.(metFieldNames{m});
        group = [];
        for f = 1:length(fieldNames)
            group = [group;metTypes.(fieldNames{f})];
        end
        if ~isempty(group)
            fraction(m,i) = sum(massFrac{i}(group));
        else
            fraction(m,i) = 0;
        end
    end
    subplot(1,length(Stoich(1,:)),i)
    pie(fraction(:,i),metFieldNames)
end

%% Growth rates
clear out
% rN = [6.5,(6.5+2.1)/2,((6.5+2.1)/2+2.1)/2, 2.1, 0.5, 0.25]/6.5*-0.2114;
model_temp = model_an;

biomassMets = Met;
biomassStoich = Stoich;
eliminateMets = {'atp_c','h2o_c'};

for i = 1:length(eliminateMets)
    elPos = find(ismember(biomassMets,eliminateMets{i}));
    biomassMets(elPos) = [];
    biomassStoich(elPos,:) = [];
end

for i = 1:length(biomassStoich(1,:))
    
    bofMetIDs = findMetIDs(model_temp,Met);
%     model_temp.lb(strmatch('EX_no3(e)',model_temp.rxns)) = rN(i);
    model_temp.S(bofMetIDs,find(model_temp.c==1)) = Stoich(:,i);
    
    out{i} = optimizeCbModel(model_temp,'max','one');
    u(i) = out{i}.f;
    
    [molarX,compMolarMass,av,Num,formula] = calculateFormula(model_temp,findMetIDs(model_temp,biomassMets),abs(biomassStoich(:,i)));
    testX = abs(sum(biomassStoich(:,i).*compMolarMass))
    
    rS = abs(out{i}.x(strmatch('EX_glc__D_e',model_temp.rxns)))*6;
    rX = abs(out{i}.f)/molarX*av(1);
    rC = abs(out{i}.x(strmatch('EX_co2_e',model_temp.rxns)))*1;
    rN = abs(out{i}.x(strmatch('EX_nh4_e',model_temp.rxns)));
    ro = 0;
    
    Ysa(i) = ro/(rS);
    Ysx(i) = rX/(rS);
    Ysc(i) = rC/(rS);
    test(i) = Ysa(i)+Ysx(i)+Ysc(i);
    
    
    Ynx(i) = rX/rN; % C-molX/molN;
    
end
R = [Ysa;Ysx;Ysc;Ynx];

%%
CSVfile = 'fluxDist_Yl';
writeFluxCsv(model_an,out,CSVfile);

%% Variations
clear deltaU deltaZ rxnID

rxnID = [1:length(model.rxns)]'*ones(1,length(u));
el = [];
for i = 2:length(Stoich(1,:))
    deltaU(:,i) = abs((out{i}.x - out{1}.x));
end
deltaU(:,1) = [];
deltaU(abs(deltaU)>100) = 0;

for i = 1:length(deltaU(:,1))
    if sum(abs(deltaU(i,:))) == 0
        el(end+1) = i;
    end
end
rxnID(el,:) = [];
deltaU(el,:) = [];
rxnID = rxnID(:,1);

deltaZ = zscore(deltaU,0,2);


%% Modification

BOF = model_an.rxns(find(model_an.c==1));
% model = changeObjective(model,BOF);

fprintf('\n \nRetrieving old and new stoichiometric coefficients\n\n')

model_analysis = model_an;

Met_ID = findMetIDs(model_analysis,Met);


fprintf('Generating optimization boundaries\n\n')
for i=1:length(Stoich(:,1))
    range_Stoech(i,:)=[min(Stoich(i,:)) max(Stoich(i,:))];
    init_Stoech(i,:)=[mean(Stoich(i,:))];
end

%%

it1 = 50; % BOF
it2 = 10; % Points

theta_0 = cell(length(Met),1);
parMetCoeff = cell(length(Met),1);
for m=1:length(Met)
    if range_Stoech(m,1) ~= range_Stoech(m,2)
        for b=1:it1 % BOF
             theta_0{m}(:,:,b)=(range_Stoech(:,1)+rand(length(range_Stoech(:,1)),1).*(range_Stoech(:,2)-range_Stoech(:,1)))*ones(1,it2);
             for p=1:it2 % Point
                 parMetCoeff{m}(p,b) = range_Stoech(m,1)+(p/it2)*(range_Stoech(m,2)-range_Stoech(m,1));
                 theta_0{m}(m,p,b) = parMetCoeff{m}(p,b);
             end
        end
    end
end
growth = cell(length(Met),1);

parfor m = 1:length(Met)
    changeCobraSolver('gurobi','all',0);
    model_analysis = model_an;
    parMet = Met(m);
    fprintf('Calculating for #%4.0f (',m)
    fprintf(char(parMet))
    fprintf(')\n\n')    
    parMetID = findMetIDs(model_analysis,parMet);
    if range_Stoech(m,1) ~= range_Stoech(m,2)
        for b =1:it1
            for p = 1:it2
                model_analysis = changeRxnMets(model_analysis,model_analysis.mets(Met_ID),model_analysis.mets(Met_ID),BOF,theta_0{m}(:,p,b));
                out=optimizeCbModel(changeObjective(model_analysis,BOF),'max','one');
                growth{m}(p,b)=out.f;
                co2_uptake{m}(p,b)=out.x(strmatch('EX_co2_e',model_analysis.rxns));
             end
        end
        for b=1:it1
           for p = 1:it2-1
               dudp{m}(p,b) = (growth{m}(p+1,b)-growth{m}(p,b))/(parMetCoeff{m}(p+1,b)-parMetCoeff{m}(p,b));
           end
        end
    else
        growth{m} = ones(it2,it1);
        parMetCoeff{m} = ones(it2,it1);
        dudp{m} = ones(it2-1,it1);
        co2_uptake{m} = ones(it2,it1);
        fprintf(char(parMet))
        fprintf(' was excluded\n\n')
    end

end
time=toc/60;
fprintf('\n\n Calculation time = %4.2f min \n\n',time)
save('BOFsensitivity_Yl.mat')