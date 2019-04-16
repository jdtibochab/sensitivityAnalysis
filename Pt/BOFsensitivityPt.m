%% Parameter sensitivity and codependence analysis
clear
load('iLB1027_lipid.mat')
model = changeObjective(iLB1027_lipid,'bof_c');
ID = 'Pt1';

%% Fixes
% M2MASN secretion
model = changeRxnBounds(model,'DM_m2masn_c',0,'b');

% Growth rate
exchange = 'EX_no3_e';
exchangeC = 'EX_co2_e';
model.lb(strmatch(exchangeC,model.rxns)) = -13.54;
[~,reqExchange] = fixGrowth(model,exchange,0.057);
model.lb(strmatch(exchange,model.rxns)) = reqExchange;

%% Read
T = readtable('PtBOF.csv');
mets = table2array(T(:,1));
Stoich = table2array(T(:,2:end));


%% Interpolation
t0 = [1	1.5	2	3	5	12	17];
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
[yieldT,yieldTc,out,~,~,~,modelsPt] = yieldAnalysis(model,mets,Stoich,'EX_co2_e');

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
T = readtable('PtBOF.csv');
mets = table2array(T(:,1));
Stoich = table2array(T(:,2:end));

%% Theoretical biomass
protFASTA = 'PtProt.fasta';
DNAFASTA = 'PtDNA.fna';
printRxnFormula(model,model.rxns(model.c==1));
[tModel,S,tS] = theoreticalBiomass(model,protFASTA,DNAFASTA,Stoich,mets,metTypes,metGroups);

tOut = optimizeCbModel(tModel);

%%
save(['BOFsensitivity_',ID])

%% Visualization
% ignoreMets = {'val__L_c','leu__L_c','gly_c','asn__L_c'};
ignoreMets = {};
costRxn = {'nadh_c -> nad_c + h_c','nadph_c -> nadp_c + h_c','ATPM'};
costMet = {'nadh_c','nadph_c','atp_c'};

sensitivityVisualization(model,mets,metTypes,metNames,growth,parMetCoeff,ID,costRxn,costMet,ignoreMets)