
met = 'Gly'
for i = 1:length(model.mets)
    a = strfind(model.metNames{i}, met);
    if a
        {model.metNames{i},model.mets{i}}
    end
end