%% Parameter sensitivity and codependence analysis
clear
load('iCHOv1.mat')
model = iCHOv1;
ID = 'CHO';
%% Read
T = readtable('CHOBOF.csv');
mets = table2array(T(:,1));
Stoich = table2array(T(:,2:end));

%% Interpolation
t0 = [0 1 2 3 4 5];
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
[yieldT,yieldTc,out,~,~,~,modelsCHO] = yieldAnalysis(model,mets,Stoich,'EX_glc__D_e');

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
protFASTA = 'CHOProt.faa';
DNAFASTA = 'CHODNA.fna';
printRxnFormula(model,model.rxns(model.c==1));
[tModel,S,tS] = theoreticalBiomass(model,protFASTA,DNAFASTA,Stoich,mets,metTypes,metGroups);

tOut = optimizeCbModel(tModel);

%%
save('BOFsensitivity_CHO.mat')

%% Visualization
% ignoreMets = {'his__L_c','phe__L_c','hista_c'};
ignoreMets = {'his__L_c'};
costRxn = {'nadh_c -> nad_c + h_c','nadph_c -> nadp_c + h_c','ATPM'};
costMet = {'nadh_c','nadph_c','atp_c'};
sensitivityVisualization(model,mets,metTypes,metNames,growth,parMetCoeff,ID,costRxn,costMet,ignoreMets)

