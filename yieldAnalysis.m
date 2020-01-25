function [T,Tc,out,fullT,fullTc,info,modelArray] = yieldAnalysis(model,mets,Stoich,substrate)
T = {};
Tc = {};
out = {};

bofMetIDs = findMetIDs(model,mets);

fullT = {};
fullTc = {};
for i = 1:length(Stoich(1,:))
    model.S(bofMetIDs,find(model.c==1)) = -Stoich(:,i);
    
    modelArray{i} = model;
    
    [tempT,tempTc,out{i},tfullT,tfullTc,info] = yieldCalculation(model,substrate);
    
    fullT(:,1) = tfullT(:,1);
    fullT(:,end+1) = tfullT(:,2);
    
    fullTc(:,1) = tfullTc(:,1);
    fullTc(:,end+1) = tfullTc(:,2);
    
end

M = fullT(:,1);
Mc = fullTc(:,1);

Y = cell2mat(fullT(:,2:end));
Yc = cell2mat(fullTc(:,2:end));

%%
M = M(find(~all(Y==0,2)));
Mc = Mc(find(~all(Yc==0,2)));

Y = Y(find(~all(Y==0,2)),:);
Yc = Yc(find(~all(Yc==0,2)),:);

T(:,1) = M;
T(:,2:length(Stoich(1,:))+1) = num2cell(Y);

Tc(:,1) = Mc;
Tc(:,2:length(Stoich(1,:))+1) = num2cell(Yc);

header = {'ID'};
for i = 1:length(Stoich(1,:))
    header{i+1} = ['t',num2str(i)];
end

T = cell2table(T,'VariableNames',header);
Tc = cell2table(Tc,'VariableNames',header);


end