function [deltaU,deltaZ,rxnID] = reactionResidualAnalysis(model,out)

rxnID = [1:length(model.rxns)]'*ones(1,length(out));
el = [];
for i = 2:length(out)
    deltaU(:,i) = abs((out{i}.x - out{1}.x));
end
deltaU(:,1) = [];
% deltaU(abs(deltaU)>100) = 0;

for i = 1:length(deltaU(:,1))
    if sum(abs(deltaU(i,:))) == 0
        el(end+1) = i;
    end
end
rxnID(el,:) = [];
deltaU(el,:) = [];
rxnID = rxnID(:,1);

deltaZ = zscore(deltaU,0,2);

end