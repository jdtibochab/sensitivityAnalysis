%%
met = 'Gly'
for i = 1:length(model.mets)
    a = strfind(model.metNames{i}, met);
    if a
        {model.metNames{i},model.mets{i}}
    end
end

%%
model = iLB1027_lipid;
model_dark.lb(4426) = 0; model_dark.ub(4426) = 1000;
model_ref = modelPA;
model = model_dark;

formulas = printRxnFormula(model_ref,'rxnAbbrList',model_ref.rxns,'printFlag',0);
solution = optimizeCbModel(model);
solution.f

for i =1:length(model_ref.rxns)
    if ~findRxnIDs(model,model_ref.rxns{i}) && all(findMetIDs(model,findMetsFromRxns(model_ref,model_ref.rxns{i})))
        model_test = addReaction(model,model_ref.rxns{i},formulas{i});
        solution = optimizeCbModel(model_test);
        solution.f
    end
    if solution.f > 0.001
        break
    end
    
end
