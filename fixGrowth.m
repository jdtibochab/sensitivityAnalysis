%%

function [model,reqExchange] = fixGrowth(model,exchange,growth)

out0 = optimizeCbModel(model);

BOFID = find(model.c==1);
BOFrxn = model.rxns(BOFID);


exID = findRxnIDs(model,exchange);


out = optimizeCbModel(model);
if out.full(exID)>0
    d = 1;
    reqExchange = model.ub(exID);
else
    d = 0;
    reqExchange = model.lb(exID);
end
i = 0;

while out.f > growth
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