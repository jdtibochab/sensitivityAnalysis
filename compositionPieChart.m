function F = compositionPieChart(model,mets,Stoich,metTypes,metGroups,ID)

N = 1;
for i = 1:length(Stoich(1,:))
    molarCoeffs = abs(Stoich(:,i));
    molarCoeffs(molarCoeffs>1) = 0;
    metIDs = findMetIDs(model,mets);
    [~,molarMass] = calculateFormula(model,metIDs,1);
    massCoeffs = molarCoeffs.*molarMass;
    totalMass(i) = sum(massCoeffs);
    massFrac{i} = massCoeffs/totalMass(i);
end

figure('visible','off','rend','painters','pos',[10 10 1500 1000])
metFieldNames = fieldnames(metGroups);
for i = 1:length(Stoich(1,:))
    for m = 1:length(metFieldNames)
        fieldNames = metGroups.(metFieldNames{m});
        group = [];
        for f = 1:length(fieldNames)
            group = [group;metTypes.(fieldNames{f})];
        end
        if ~isempty(group)
            fraction(m,i) = sum(massFrac{i}(group));
        else
            fraction(m,i) = 0;
        end
    end
    subplot(N,length(Stoich(1,:))/N,i)
    pie(fraction(:,i),metFieldNames)
end
saveas(gcf,['compositionPieChart_',ID],'svg')

%%

F(:,1) = metFieldNames;
F(:,2:length(fraction(1,:))+1) = num2cell(fraction);

header = {'Group'};
for i = 1:length(fraction(1,:))
    header{i+1} = ['t',num2str(i)];
end

F = cell2table(F,'VariableNames',header);

end