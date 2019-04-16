clear
% Pt
tic
fprintf('\n\n Calculating for CHO TEST1\n\n')

load('iCHOv1.mat')
model_an = iCHOv1;

global Met_ID

%%
% T=readtable('BOFcho.xls');
% mets = T.Var1;
% metIDs = findMetIDs(model_an,mets);
% Stoich = [T.x0,T.x1,T.x2,T.x3,T.x4,T.x5];
% [~,molarMass] = calculateFormula(model_an,metIDs,1);
% 
% for i=1:length(Stoich(1,:))
%     massCoeff = Stoich(:,i);
%     total = sum(massCoeff);
%     newMassCoeff = massCoeff/total;
%     newStoichCoeff = newMassCoeff./molarMass;
%     Stoich(:,i) = newStoichCoeff;
% end
% 
% [model_an,prodCapMet,notProdCapMet,Stoich2,theta_0] = metTest(model_an,mets,Stoich);
% save('CHOdata.mat')


load ('CHOdata.mat')
Stoich = [Stoich2;ones(1,length(Stoich(1,:)))*28.93;ones(1,length(Stoich(1,:)))*24.7485];
Met = [prodCapMet;'atp_c';'h2o_c'];
Stoich = Stoich*-1;

%%
colors = {'y','m','c','r','g','b'};
t0 = [0 1 2 3 4 5];
% t = [0.5 1 1.5 2 3.5 5];
t = [min(t0):(max(t0)-min(t0))/5:max(t0)];
[intStoich,rsq] = interpolateStoich(Stoich,t0,t);
Stoich = intStoich;
%% Biomass composition change
N = 1;
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

%% Modification
clear out
biomassReactionIDs = [strmatch('BIOMASS_',model_an.rxns);strmatch('bof',model_an.rxns)];
biomassReactions = model_an.rxns(biomassReactionIDs);

model_an.ub(biomassReactionIDs(2)) = 0;

BOF = model_an.rxns(find(model_an.c==1));
% model = changeObjective(model,BOF);

CSVfile = 'fluxDist_CHO';
model_temp = model_an;
for i = 1:length(Stoich(1,:))
    bofMetIDs = findMetIDs(model_temp,Met);
%     model_temp.lb(strmatch('EX_no3_e',model_temp.rxns)) = rN(i);
    model_temp.S(bofMetIDs,find(model_temp.c==1)) = Stoich(:,i);
    modelsPt{i} = model_temp;
    out{i} = optimizeCbModel(model_temp,'max','one');
    u(i) = out{i}.f;
end

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
fprintf('\n \nRetrieving old and new stoichiometric coefficients\n\n')

model_analysis = model_an;

Met_ID=findMetIDs(model_analysis,Met);

fprintf('Generating optimization boundaries\n\n')
for i=1:length(Stoich(:,1))
    range_Stoech(i,:)=[min(Stoich(i,:)) max(Stoich(i,:))];
    init_Stoech(i,:)=[mean(Stoich(i,:))];
end

%%
clear out
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

% parpool(12)
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
save('BOFsensitivity_CHO.mat')