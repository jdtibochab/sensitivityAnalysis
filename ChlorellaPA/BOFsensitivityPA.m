
%% Parameter sensitivity and codependence analysis
clear
load('ModelPA')
model = ModelPA{1};
ID = 'PA';

%% Fixes
model = starchCorr(model,'starch300_h','DM_o2D(u)',0.0355);

%% Read
T = readtable('chlorellaPABOF.csv');
mets = table2array(T(:,1));
Stoich = table2array(T(:,2:end));

%% Interpolation
t0 = [3 4 5 6 7 8];
t = [min(t0):(max(t0)-min(t0))/5:max(t0)];
oldStoich = Stoich(:,1:6);
newStoich = Stoich(:,7:end);
[oldStoich,rsq] = interpolateStoich(oldStoich,t0,t);
[newStoich,rsq] = interpolateStoich(newStoich,t0,t);
% Stoich = [oldStoich,newStoich]; 
Stoich = (oldStoich+newStoich)/2;

% Growth interpolation
uE = [0.0274847619	0.0279375039	0.0227957421	0.0207400832	0.0189629058	0.0175230806];
uE = interpolateStoich(uE,t0,t);

%% Test for production capacity
[model,mets,notProdCapMet,Stoich,~] = metTest(model,mets,Stoich);

%% Metabolite classification
[metTypes,metGroups,metNames] = findMetType(model,mets);

%% Biomass composition change
F = compositionPieChart(model,mets,Stoich,metTypes,metGroups,ID);

%% Yield changes
[a,b] = exchangeSingleModel(model);
[yieldT,yieldTc,out,~,~,~,modelsPA] = yieldAnalysis(model,mets,Stoich,'EX_co2(e)');

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
protFASTA = 'CvProt.fasta';
DNAFASTA = 'CvDNA.fasta';
printRxnFormula(model,model.rxns(model.c==1));
[tModel,S,tS] = theoreticalBiomass(model,protFASTA,DNAFASTA,Stoich,mets,metTypes,metGroups);

tOut = optimizeCbModel(tModel);

%%
save(['BOFsensitivity_',ID])

%% Visualization
% ignoreMets = {'starch','dgtp_c','gtp_c','datp_c'};
ignoreMets = {'starch'};
costRxn = {'nadh_c -> nad_c + h_c','nadph_c -> nadp_c + h_c','ATPM(NGAM)'};
costMet = {'nadh_c','nadph_c','atp_c'};
sensitivityVisualization(model,mets,metTypes,metNames,growth,parMetCoeff,ID,costRxn,costMet,ignoreMets)


