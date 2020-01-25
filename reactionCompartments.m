
function [rxnComps,allComps,comps] = reactionCompartments(model)

for i=1:length(model.mets)
    met = model.mets{i};
    if strcmp(met(end),']');
        allComps{i,1} = met(end-1);
    else
        allComps{i,1} = met(end);
    end
    
end

% Correct extra compartments
incorrectComps = {'o','w','s','p','i','a'};
correctComps = {'n','v','r','m','g','m'};
for i = 1:length(allComps)
    compPos = strmatch(allComps{i},incorrectComps);
    if ~isempty(compPos)
        allComps(i) = correctComps(compPos);
    end
end

comps = unique(allComps);

rxnComps = zeros(length(model.rxns),1);
for r = 1:length(model.rxns)
    metsInRxn = find(model.S(:,r)~=0);
    compsInRxn = allComps{metsInRxn};
    numComps = length(unique(compsInRxn));
    if numComps == 1
        comp = unique(compsInRxn);
        compID = find(ismember(comps,comp));
        rxnComps(r) = compID;
    end

end


% cellstr(cellfun(@(v)v(end),model.mets)))