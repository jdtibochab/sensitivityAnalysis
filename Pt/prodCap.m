
%% Test production capability one by one

function [model,prodCapMets,notProdCapMets,Stoech] = prodCap(model,Met,Stoech)
% Pt
biomassReactionIDs = [strmatch('BIOMASS_',model.rxns);strmatch('bof',model.rxns)];
model_an.lb(biomassReactionIDs(1:end-1)) = 0;
model_an.ub(biomassReactionIDs(1:end-1)) = 0;

BOFid = find(model.c==1);

fprintf('\n\n Calculating for Pt\n\n')

% load('ptData.mat')
metIDs = find(ismember(model.mets,Met));

size(metIDs)
size(Met)

prodCapMetIDs = [];
prodCapMetNum = [];
prodCapMets = {};
notProdCapMetIDs = [];
notProdCapMetNum = [];
notProdCapMets = {};
excMetIDs = [];
excMetNum = [];

%% Modification

biomassReactionIDs = [strmatch('BIOMASS_',model.rxns);strmatch('bof',model.rxns)];
modelmets = 0;
for i = 1:length(biomassReactionIDs)
    k=length(modelmets);
    modelmets_temp = find(model.S(:,biomassReactionIDs(i))~=0);
    modelmets(k:k-1+length(modelmets_temp)) = modelmets_temp;
end
modelmets = sort(modelmets');
k=0;
e=0;
for m = 1:length(Met)
    if sum(ismember(model.mets(modelmets),Met(m))) == 0
        k = k+1;
        addMetIDs(k,1) = metIDs(m);
        addMetNum(k,1) = m;
    else
        e = e+1;
        excMetIDs(e,1) = metIDs(m);
        excMetNum(e,1) = m;
        fprintf('\n')
        fprintf(char(Met(m)))
        fprintf(' was excluded \n')
    end
end
p=0;
n=0;

metCoeff = -mean(Stoech')';

for j = 1:length(addMetIDs)
    fprintf('\n')
    fprintf(char(model.mets(addMetIDs(j))))
    temp_model = model;
    temp_model.S(find(temp_model.S(:,BOFid)<0),BOFid) = temp_model.S(find(temp_model.S(:,BOFid)<0),BOFid)/(1+metCoeff(addMetNum(j)));
    temp_model.S(addMetIDs(j),BOFid) = metCoeff(addMetNum(j));

    out = optimizeCbModel(temp_model,'max');
    if out.stat == 1 && out.f > 1e-4
       p = p+1;
       prodCapMetIDs(p,1) = addMetIDs(j);
       prodCapMetNum(p,1) = addMetNum(j);
       fprintf(' yes\n')
    else
       n=n+1;
       notProdCapMetIDs(n,1) = addMetIDs(j);
       notProdCapMetNum(n,1) = addMetNum(j);
       fprintf(' no\n')
       notProdCapMets(n,1) = Met(notProdCapMetNum(n,1));
       fprintf(char(Met(notProdCapMetNum(n,1))));
       fprintf(' was discarded \n')
    end
end


if n
    for i = 1 : length(notProdCapMetNum)
        j = length(notProdCapMetNum)-i+1;
        Stoech(notProdCapMetNum(j),:) = [];
    end
end

prodCapMetNum = sort([excMetNum;prodCapMetNum]);

prodCapMets = Met(prodCapMetNum);

size(metCoeff)
prodCapMetNum

model.S(find(ismember(model.mets,prodCapMets)),BOFid) = metCoeff(prodCapMetNum)

end
