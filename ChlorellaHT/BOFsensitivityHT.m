%% Parameter sensitivity and codependence analysis
clear
load('miCZ_H')
ID = 'HT';

%%
model = modelDa;

%% Fixes
model = starchCorr(model,'starch300_h','EX_no3(e)',0.0355);
% model = addReaction(model,'ICITht','reactionFormula','icit_h <=> icit_c');
% model = addReaction(model,'ICITht','reactionFormula','icit_h  -> glx_h + succ_h');
% model = addReaction(model,'ICITht','reactionFormula','mal_L_c + icit_h  <=> mal_L_h + icit_c');
model = addReaction(model,'ICITht','reactionFormula','akg_c + icit_h  <=> akg_h + icit_c','subSystem','Transport');
% model = addReaction(model,'ICITht','reactionFormula','icit_h + oaa_c  <=> icit_c + oaa_h');

model = changeRxnBounds(model,'DM_icit(h)',0,'b');

%% Read composition data
T = readtable('chlorellaHTBOF.csv');
mets = table2array(T(:,1));
Stoich = table2array(T(:,2:end));

%% Interpolation
t0 = [4.24	4.75	5.24	5.75	6.25];
t = [min(t0):(max(t0)-min(t0))/5:max(t0)];
oldStoich = Stoich(:,1:5);
newStoich = Stoich(:,6:end);
[oldStoich,rsq] = interpolateStoich(oldStoich,t0,t);
[newStoich,rsq] = interpolateStoich(newStoich,t0,t);
% Stoich = [oldStoich,newStoich]; 
Stoich = (oldStoich+newStoich)/2;

% Growth interpolation
uE = [0.0256	0.0263	0.0254	0.0243	0.0236];
uE = interpolateStoich(uE,t0,t);

%% Test for production capacity
[model,mets,notProdCapMet,Stoich,~] = metTest(model,mets,Stoich);

%% Metabolite classification
[metTypes,metGroups,metNames] = findMetType(model,mets);

%% Biomass composition change
F = compositionPieChart(model,mets,Stoich,metTypes,metGroups,ID);

%% Yield changes
[a,b] = exchangeSingleModel(model);
[yieldT,yieldTc,out,~,~,~,modelsHT] = yieldAnalysis(model,mets,Stoich,'EX_glc-A(e)');

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
ignoreMets = {'starch','dgtp_c','gtp_c','datp_c','chla_u','chlb_u'};
costRxn = {'nadh_c -> nad_c + h_c','nadph_c -> nadp_c + h_c','ATPM(NGAM)'};
costMet = {'nadh_c','nadph_c','atp_c'};

sensitivityVisualization(model,mets,metTypes,metNames,growth,parMetCoeff,ID,costRxn,costMet,ignoreMets)

