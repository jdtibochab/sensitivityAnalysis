load('miCZ_H.mat');

model = modelDa;


metIDs = findMetIDs(model,mets);

incorrectMets = mets(metIDs==0);

for m = 1:length(incorrectMets)
    met = incorrectMets{m};
    str = strtok(incorrectMets{m},'_');

    metID = strmatch(str,model.mets);
    
    for i = 1:length(metID)
        incMet = model.mets{metID(i)};
        [temp,b] = strtok(incMet,'-');
        if length(temp) == 3
            [a,~] = strtok(b,'-');
            corrMet = [temp,'__',a];
        else
            corrMet = incMet;
        end
        model.mets{metID(i)} = corrMet;
    end
    
    
end