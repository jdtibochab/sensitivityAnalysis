%% Parameter sensitivity and codependence analysis
clear
tic

% HETEROTROPHY
load('ModelHTold')
Model = ModelH;
load('ModelHT')

load('modelHT')

fprintf('\n\n Calculating for heterotrophic growth\n\n')

model_an = [Model,ModelH];
for m = 1:length(model_an)
   model_an{m} = starchCorr(model_an{m},'starch300_h','EX_no3(e)',0.02);
%    model_an{m} = changeRxnBounds(model_an{m},'EX_starch(h)',0,'l');
%    model_an{m} = changeRxnBounds(model_an{m},'EX_co2(e)',0,'l');
%    model_an{m} = changeRxnBounds(model_an{m},'DM_icit(h)',0,'b');
end


solutionH = cell(1,length(model_an));
growth_initial = zeros(length(model_an));
for j=1:length(model_an)
    model=model_an{1,j};
    printRxnFormula(model,model.rxns(find(model.c==1)));
    solutionH{j}=optimizeCbModel(model);
    growth_initial(j) = solutionH{j}.f;
end
    

fprintf('\n \nRetrieving old and new stoichiometric coefficients\n\n')

Srxn2 = cell(1,length(model_an));
Met2 = cell(1,length(model_an));
for j=1:length(model_an)
    model=model_an{1,j};
    Srxn = full(model.S(:,find(model.c==1)));
        Srxn2{j} = Srxn(find(Srxn~=0));
    Met=model.mets(find(Srxn~=0));
        Met2{j} = Met;
end

lengths = zeros(length(model_an));
for i=1:length(model_an)
    lengths(i)=length(Srxn2{i});
end
maxLength = max(lengths);

%% Vector length correction
for i=1:6
    for j=1:maxLength
       a=strcmp(Met2{i}(j),Met2{7}(j));
       if a == 0
           Met2{i}=[Met2{i}(1:j-1);Met2{7}(j);Met2{i}(j:end)];
           Srxn2{i}=[Srxn2{i}(1:j-1);0;Srxn2{i}(j:end)];
       end
    end   
end

Stoich = cell2mat(Srxn2);
%%
t0 = [4.24	4.75	5.24	5.75	6.25];
% t = [4.24 4.5 5 5.5 6 6.5];
t = [min(t0):(max(t0)-min(t0))/5:max(t0)];
oldStoich = Stoich(:,1:5);
newStoich = Stoich(:,6:end);
[oldStoich,rsq] = interpolateStoich(oldStoich,t0,t);
[newStoich,rsq] = interpolateStoich(newStoich,t0,t);
% Stoich = [oldStoich,newStoich]; 
% Stoich = (oldStoich+newStoich)/2;

%% Biomass composition change
N = 1;
mets = Met2{1};
for i = 1:length(Stoich(1,:))
    molarCoeffs = abs(Stoich(:,i));
    molarCoeffs(molarCoeffs>1) = 0;
    metIDs = findMetIDs(model_an{1},mets);
    [~,molarMass] = calculateFormula(model,metIDs,1);
    massCoeffs = molarCoeffs.*molarMass;
    totalMass(i) = sum(massCoeffs);
    massFrac{i} = massCoeffs/totalMass(i);
end
[metTypes,mets] = findMetType(model,mets);

metGroups.CB = {'CB'};
metGroups.rest = {'rest'};
metGroups.Nuc = {'Nuc'};
metGroups.AA = {'AA'};
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
    subplot(N,length(Stoich(1,:))/N,i)
    pie(fraction(:,i),metFieldNames)
end
%%
model_analysis = modelHT;

Met_ID=findMetIDs(model_analysis,Met);

fprintf('Generating optimization boundaries\n\n')

range_Stoech = zeros(length(Met),2);
for i=1:length(Met)
    range_Stoech(i,:)=[min(Stoich(i,:)) max(Stoich(i,:))];
end

%%
clear out

model_temp = model_an{1};

biomassMets = Met;
biomassStoich = Stoich;
eliminateMets = {'h_c','pi_c','atp_c','h2o_c','adp_c'};

for i = 1:length(eliminateMets)
    elPos = find(ismember(biomassMets,eliminateMets{i}));
    biomassMets(elPos) = [];
    biomassStoich(elPos,:) = [];
end

CSVfile = 'fluxDist_PA';
for i = 1:length(biomassStoich(1,:))
    bofMetIDs = findMetIDs(model_temp,biomassMets);
%     model_temp.lb(strmatch('EX_no3(e)',model_temp.rxns)) = rN(i);
    model_temp.S(bofMetIDs,find(model_temp.c==1)) = biomassStoich(:,i);
    modelsPA{i} = model_temp;
    out{i} = optimizeCbModel(model_temp,'max','one');
    u(i) = out{i}.f;
    
    [molarX,~,av] = calculateFormula(model_temp,findMetIDs(model_temp,biomassMets),abs(biomassStoich(:,i)));
    
    rS = abs(out{i}.x(strmatch('EX_glc-A(e)',model_temp.rxns)))*6;
    rX = abs(out{i}.f)/molarX*av(1);
    rC = abs(out{i}.x(strmatch('EX_co2(e)',model_temp.rxns)))*1;
    rN = abs(out{i}.x(strmatch('EX_no3(e)',model_temp.rxns)));
    ro = abs(out{i}.x(strmatch('DM_icit(h)',model_temp.rxns)))*6;
    
    Ysa(i) = ro/(rS);
    Ysx(i) = rX/(rS);
    Ysc(i) = rC/(rS);
    test(i) = Ysa(i)+Ysx(i)+Ysc(i);
    
    Ynx(i) = rX/rN; % C-molX/molN;
    
end
R = [Ysa;Ysx;Ysc;Ynx];

%%
CSVfile = 'fluxDist_HT';
writeFluxCsv(model_temp,out,CSVfile);

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
%%
it1 = 100; % BOF
it2 = 10; % Points

theta_0 = cell(length(Met),1);
parMetCoeff = cell(length(Met),1);
for k=1:length(Met)
    if range_Stoech(k,1) ~= range_Stoech(k,2)
        for j=1:it1 % BOF
             theta_0{k}(:,:,j)=(range_Stoech(:,1)+rand(length(range_Stoech(:,1)),1).*(range_Stoech(:,2)-range_Stoech(:,1)))*ones(1,it2);
             for i=1:it2 % Point
                 parMetCoeff{k}(i,j) = range_Stoech(k,1)+(i/it2)*(range_Stoech(k,2)-range_Stoech(k,1));
                 theta_0{k}(k,i,j) = parMetCoeff{k}(i,j);
             end
        end
    end
end
growth = cell(length(Met),1);

parfor k=1:length(Met)
    changeCobraSolver('gurobi','all',0);
    parMet = Met(k);
    fprintf('Calculating for #%4.0f (',k)
    fprintf(char(parMet))
    fprintf(')\n\n')    
    parMetID = findMetIDs(model_analysis,parMet);
    if range_Stoech(k,1) ~= range_Stoech(k,2)
        for j =1:it1
            for i = 1:it2
                model_analysis = changeRxnMets(model_analysis,model_analysis.mets(Met_ID),model_analysis.mets(Met_ID),'Biomass_Cvu_hetero-',theta_0{k}(:,i,j));
                out=optimizeCbModel(changeObjective(model_analysis,'Biomass_Cvu_hetero-'),'max','one');
                growth{k}(i,j)=out.f;
                co2_uptake{k}(i,j)=out.x(strmatch('EX_glc-A(e)',model_analysis.rxns));
             end
        end
        for j=1:it1
           for i = 1:it2-1
               dudp{k}(i,j) = (growth{k}(i+1,j)-growth{k}(i,j))/(parMetCoeff{k}(i+1,j)-parMetCoeff{k}(i,j));
           end
        end
    else
        growth{k} = ones(it2,it1);
        parMetCoeff{k} = ones(it2,it1);
        dudp{k} = ones(it2-1,it1);
        co2_uptake{k} = ones(it2,it1);
        fprintf(char(parMet))
        fprintf(' was excluded\n\n')
    end

end
time=toc/60;
fprintf('\n\n Calculation time = %4.2f min \n\n',time)
save('BOFsensitivity_HT2.mat')