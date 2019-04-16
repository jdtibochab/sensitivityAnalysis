%% Parameter sensitivity and codependence analysis
clear
load('iMM904_genes4.19.2016.mat')
model = modelSc2;
ID = 'Sc';

%% Fixes
model = changeRxnBounds(model,'EX_nh4_e',-1000,'l');
optimizeCbModel(model)

%% Read
T = readtable('ScBOF.csv');
mets = table2array(T(:,1));
Stoich = table2array(T(:,2:end));

%% Interpolation
t0 = [8 12 16]; % h
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
[yieldT,yieldTc,out,~,~,~,modelsSc] = yieldAnalysis(model,mets,Stoich,'EX_glc__D_e');

%% Flux CSV files
CSVfile = ['fluxDist_',ID];
writeFluxCsv(model,out,CSVfile);

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
protFASTA = 'ScProt.faa';
DNAFASTA = 'ScDNA.fna';
printRxnFormula(model,model.rxns(model.c==1));
[tModel,S,tS] = theoreticalBiomass(model,protFASTA,DNAFASTA,Stoich,mets,metTypes,metGroups);

tOut = optimizeCbModel(tModel);

%%
save(['BOFsensitivity_',ID])

%% Visualization
ignoreMets = {};
costRxn = {'nadh_c -> nad_c + h_c','nadph_c -> nadp_c + h_c','ATPM'};
costMet = {'nadh_c','nadph_c','atp_c'};

sensitivityVisualization(model,mets,metTypes,metNames,growth,parMetCoeff,ID,costRxn,costMet,ignoreMets)