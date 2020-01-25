%% Parameter sensitivity and codependence analysis
clear
load('iYL_ready.mat')
model = changeObjective(modelD,'xBIOMASS');
ID = 'Yl';

%% Read
T = readtable('YlBOF.csv');
mets = table2array(T(:,1));
Stoich = table2array(T(:,2:end));

%% Interpolation
t0 = [0.5 1.625 3.62];
t = [min(t0):(max(t0)-min(t0))/5:max(t0)];
[intStoich,rsq] = interpolateStoich(Stoich,t0,t);

Stoich = intStoich;

%% Test for production capacity
[model,mets,notProdCapMet,Stoich,~] = metTest(model,mets,Stoich);

%% Metabolite classification
[metTypes,metGroups,metNames] = findMetType(model,mets);

%% Biomass composition change
F = compositionPieChart(model,mets,Stoich,metTypes,metGroups,ID);

%% Yield changes
[a,b] = exchangeSingleModel(model);
[yieldT,yieldTc,out,~,~,~,modelsYl] = yieldAnalysis(model,mets,Stoich,'EX_glc__D_e');

%% Nucleus
Nmets = {};

for m = 1:length(model.mets)
    met = model.mets{m};
    if strcmp('n',met(end))
        Nmets{end+1,1} = met;
    end
end

NmetIDs = findMetIDs(model,Nmets);

Nrxns = {};
for m = 1:length(Nmets)
    NrxnID = find(model.S(NmetIDs(m),:));
    for r = 1:length(NrxnID)
        if out{3}.x(NrxnID(r))
            Nrxns{end+1} = Nrxn;
        end
    end
end

%% Flux CSV files
load('allModels.mat'); [compModel,H,NH] = homogenizeModel(model,modelSc,unique(model.subSystems));
CSVfile = ['fluxDist_',ID];
writeFluxCsv(compModel,out,CSVfile);

%% Variations
[deltaU,deltaZ,rxnID] = reactionResidualAnalysis(model,out);

%% BOF sensitivity
itBOF = 50;
itPoint = 10;

calcMets = mets;
% calcMets = mets(metTypes.AA);

tic;
[growth,parMetCoeff,dudp] = BOFsensitivity(model,mets,calcMets,Stoich,itBOF,itPoint,12);
time = toc/60;

%% Theoretical biomass
protFASTA = 'YlProt.fasta';
DNAFASTA = 'YlDNA.fna';
printRxnFormula(model,model.rxns(model.c==1));
[tModel,S,tS] = theoreticalBiomass(model,protFASTA,DNAFASTA,Stoich,mets,metTypes,metGroups);

tOut = optimizeCbModel(tModel);

%% 
save(['BOFsensitivity_',ID])

%% Visualization
ignoreMets = {'pg_p'};
% ignoreMets = {};
costRxn = {'nadh_c -> nad_c + h_c','nadph_c -> nadp_c + h_c','227'};
costMet = {'nadh_c','nadph_c','atp_c'};

sensitivityVisualization(model,mets,metTypes,metNames,growth,parMetCoeff,ID,costRxn,costMet,ignoreMets)
