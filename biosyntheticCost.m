function cost = biosyntheticCost(model,mets,costRxn)
%%
metIDs = findMetIDs(model,mets);

bof = model.rxns(model.c==1);
bofID = findRxnIDs(model,bof);

costRxnID = findRxnIDs(model,costRxn);
if ~costRxnID
    model = addReaction(model,'COST','reactionFormula',costRxn,'printLevel',0);
    costRxn = model.rxns(end);
end

%%
out0 = optimizeCbModel(model);
maxGrowthRate = out0.f;
fixedGrowthRate = 0.95 * maxGrowthRate;

%%
model = changeRxnBounds(model,bof,fixedGrowthRate,'b');
model = changeRxnBounds(model,costRxn,0,'l');
model = changeRxnBounds(model,costRxn,100000,'u');

model = changeObjective(model,costRxn);

out = optimizeCbModel(model);
costCons = out.f;

%%
for m = 1:length(mets)
    modelDM = addDemandReaction(model,mets{m},0);
    DMrxn = ['DM_',char(mets{m})];
    DMflux = 0.001;
    modelDM = changeRxnBounds(modelDM,DMrxn,DMflux,'b');
    outDM = optimizeCbModel(modelDM);
    costConsDM = outDM.f;
    
    cost(m) = (costCons - costConsDM)/DMflux; % mmol[Cost Met]/mmolMet
end

end