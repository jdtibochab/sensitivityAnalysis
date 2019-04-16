function C=calculateConectivity(model,mets)
metIDs = findMetIDs(model,mets);

% out = optimizeCbModel(model);
% rxnIDs = find(out.full~=0); % Only those reactions carrying flux

rxnIDs = findRxnIDs(model,model.rxns);

% 
% length(mets)
% length(metIDs)
for i=1:length(mets)
   A = find(model.S(metIDs(i),rxnIDs)~=0);
   C(i,1) = length(A);
end

end