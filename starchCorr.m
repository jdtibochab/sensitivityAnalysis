%% Starch correction
function model = starchCorr(model,starchMet,ex,target)

out0 = optimizeCbModel(model);

starchMetID = findMetIDs(model,starchMet);

BOFID = find(model.c==1);
BOFrxn = model.rxns(BOFID);

stoechCoeff = model.S(starchMetID,BOFID);
newStoechCoeff = stoechCoeff*4/300;

model.S(starchMetID,BOFID) = newStoechCoeff;

exID = findRxnIDs(model,ex);


out = optimizeCbModel(model);
if out.full(exID)>0
    d = 1;
    reqExchange = model.ub(exID);
else
    d = 0;
    reqExchange = model.lb(exID);
end
i = 0;
while out.f > target
    i = i+1;
    reqExchange = reqExchange*0.95;
    if d
        model.ub(exID) = reqExchange;
    else
        model.lb(exID) = reqExchange;
    end
    out = optimizeCbModel(model);
end

end