function [T,Tc,out,fullT,fullTc,info] = yieldCalculation(model,substrate)

bofID = find(model.c==1);
metIDs = find(model.S(:,bofID)<0);
mets = model.mets(metIDs);
Stoich = model.S(metIDs,bofID);

notBiomassMet = {'atp_c','adp_c','h2o_c','h_c','pi_c'};

notBiomassPos = find(ismember(mets,notBiomassMet));

biomassMet = mets;
biomassMet(notBiomassPos) = [];
biomassStoich = Stoich;
biomassStoich(notBiomassPos,:) = [];
biomassMetIDs = findMetIDs(model,biomassMet);

Stoich = biomassStoich;
mets = biomassMet;
metIDs = biomassMetIDs;

%%
exchangeRxnIDs = find(findExcRxns(model));
exchangeRxns = model.rxns(exchangeRxnIDs);

notExchangeMet = {'biomass'};

for e = 1:length(exchangeRxnIDs)
    exchangeMets(e,1) = findMetsFromRxns(model,exchangeRxns{e});
    exchangeMetIDs(e,1) = findMetIDs(model,exchangeMets{e});
end

notExchangePos = find(ismember(exchangeMets,notExchangeMet));
exchangeMets(notExchangePos) = [];
exchangeMetIDs(notExchangePos) = [];
exchangeRxns(notExchangePos) = [];
exchangeRxnIDs(notExchangePos) = [];

metFormulas = model.metFormulas(exchangeMetIDs);
for e = 1:length(exchangeRxnIDs)
    convFac(e,1) = numAtomsOfElementInFormula(char(metFormulas(e)),'C');
end

carbonYields = find(convFac>0);
convFac(convFac==0) = 1;

%%
substrateRxnID = findRxnIDs(model,substrate);
substrateMetID = find(model.S(:,substrateRxnID)~=0);
convFacSubs = numAtomsOfElementInFormula(char(model.metFormulas(substrateMetID)),'C');

%%
out = optimizeCbModel(model,'max','one');

[molarX,compMolarMass,av,~,biomassFormula] = calculateFormula(model,metIDs,abs(Stoich));
testX = abs(sum(Stoich.*compMolarMass));

if abs(testX-1) > 0.01
    warning('Mass balance not met')
end

rS = abs(out.x(strmatch(substrate,model.rxns))*convFacSubs);
rX = out.f/molarX*av(1);
Ysx = rX/rS;

for j = 1:length(exchangeMets)
    r(j,1) = out.x(exchangeRxnIDs(j))*convFac(j);
    Y(j,1) = r(j,1)/rS;
end

%%

convFac = convFac(carbonYields);
convFac = [av(1);convFac(~all(Y(carbonYields,:)==0,2))];

metFormulas = metFormulas(carbonYields);
formula = [biomassFormula;metFormulas(~all(Y(carbonYields,:)==0,2))];

Mc = ['X';exchangeMets(carbonYields)];
Yc = [Ysx;Y(carbonYields,:)];

fullTc(:,1) = Mc;
fullTc(:,2:length(Stoich(1,:))+1) = num2cell(Yc);

M = ['X';exchangeMets];
Y = [Ysx;Y];
fullT(:,1) = M;
fullT(:,2:length(Stoich(1,:))+1) = num2cell(Y);

M = M(find(~all(Y==0,2)));
Mc = Mc(find(~all(Yc==0,2)));

Y = Y(find(~all(Y==0,2)),:);
Yc = Yc(find(~all(Yc==0,2)),:);

T(:,1) = M;
T(:,2:length(Stoich(1,:))+1) = num2cell(Y);

Tc(:,1) = Mc;
Tc(:,2:length(Stoich(1,:))+1) = num2cell(Yc);

%%
info = Tc(:,1);
info(:,2) = num2cell(convFac);
info(:,3) = formula;
end