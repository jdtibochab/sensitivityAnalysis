%% Clustergrams
clear
loadFileID = 'BOFsensitivity_';

models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};
% models = {'PA'};

DELTA = [];
for K = 1:length(models)
    markerSubs = [];
    NuMS = [];
    deltaZ = {};
    rxnID = [];
    clear model
    
    load([loadFileID,models{K}])
    Zscore{K} = deltaZ;
    modelStructs{K} = model;
    reactionIDs{K} = rxnID;
    FBA.(models{K}) = out;
    
    for r = 1:length(rxnID)
        metabolites{K} = findMetsFromRxns(modelStructs{K},rxnID);
    end
    
    DELTA = [DELTA;Zscore{K}];
    
    colLab = {'1','2','3','4','5'};
    cluster = clustergram(Zscore{K},'Colormap',redbluecmap,'ColumnLabels',colLab);
    modelClustergrams.(models{K}) = cluster;
    % Colormap
    [originalG_subs{K},G_subs{K}] = generalcompartmentsModel(model.subSystems(rxnID));
    
    uG_subs{K} = unique(G_subs{K});
    
    for i = 1:length(G_subs{K})
        markerSubs(i,1) = find(ismember(uG_subs{K},G_subs{K}{i}));
    end
    uMarkerSubs = unique(markerSubs);
    
%     
    % Eliminate minor subsystems
    el = [];
    for i = 1:length(uMarkerSubs)
        NuMS(i) = length(find(markerSubs == uMarkerSubs(i)));
    end
    weight = NuMS/sum(NuMS);
    elSubs = find(weight<0.01);
    
    elGsubs = [];
    for i = 1:length(elSubs)
        elGsubs = [elGsubs;find(ismember(G_subs{K},uG_subs{K}(elSubs(i))))];
    end
    
    uG_subs{K}(elSubs) = [];
    G_subs{K}(elGsubs) = [];
end
close all hidden
save('clustergrams.mat')

%% Subsystem heatmaps
clear colors
% colorIDs = {'SaddleBrown';'LightYellow';'Gold';'Khaki';'Plum';'LightGreen';'White'};
% for c = 1:length(colorIDs)
%     colors(c,:) = rgb(colorIDs{c});
% end

colorIDs = {'#ff7f00';'#e41a1c';'#984ea3';'#ffffb3';'#377eb8';'#4daf4a';'#a65628'};
colors = hex2rgb(colorIDs);

allSubsystems = [];
for K = 1:length(models)
    allSubsystems=[allSubsystems;uG_subs{K}];
end
allSubsystems = unique(allSubsystems);

eliminateSubsystems = {'Other','Cell envelope','Biomass, maintenance, demand, and expectral decomposition','Transport reactions','Terpenoid biosynthesis'};

allSubsystems(find(ismember(allSubsystems,eliminateSubsystems))) = [];
allSubsystems{end+1} = 'Other';

% allSubsystems = {'Amino acid metabolism','Carbohydrate metabolism','Cell envelope','Energy metabolism',...
%     'Fatty acid metabolism','Metabolis of cofactors and vitamins','Metabolism of other amino acids','Nucleotide metabolism','Other'};

figure
allMarkers = [];
for K = 1:length(models)
    newMarkerSubs{K} = [];
    uMarkerSubs = [];
    for i = 1:length(G_subs{K})
        index = find(ismember(allSubsystems,G_subs{K}{i}));
        if ~isempty(index)
            newMarkerSubs{K}(i,1) = index;
        else
            newMarkerSubs{K}(i,1) = length(allSubsystems);
        end
    end
    uMarkerSubs = unique(newMarkerSubs{K});
    allMarkers = [allMarkers;newMarkerSubs{K}];
%     subplot(2,4,K)
%     modelHM.(models{K}) = HeatMap([newMarkerSubs,newMarkerSubs],'Colormap',colormap(colorcube(length(allSubsystems))),'Symmetric',false,'DisplayRange',length(allSubsystems));
%     modelHM.(models{K}) = heatmap([newMarkerSubs,newMarkerSubs]);
%     caxis([1 length(allSubsystems)])
%     colorbar
end

for K = 1:length(models)
    figure
    modelHM.(models{K}) = heatmap([newMarkerSubs{K},newMarkerSubs{K}]);
    colormap(colors)
    caxis([1 length(allSubsystems)])
    saveas(gcf,['subsystemsHM_',models{K}],'svg')
end

%% Check other subs
for K = 1:length(models)
    otherSubs{K} = {};
    for i = 1:length(G_subs{K})
        if strcmp('Other',G_subs{K}{i})
            otherSubs{K}{end+1} = originalG_subs{K}{i};
        end
    end
end

%% Clustergrams
for K = 1:length(models)
    modelClustergrams.(models{K}).plot;
%     saveas(gcf,['heatmap',models{K}],'png')
end

%% Histogram
figure
counts = hist(allMarkers,length(allSubsystems));
[counts,index] = sort(counts);
pie(counts)

colormap(colors(index,:))

%% Compartments load (pie charts)
colorIDs = {'Black';'LightBlue';'Silver';'Green';'DarkSlateGray';'Pink';'Tomato';'Brown';'PaleGreen';'SaddleBrown';'GoldenRod';'White'};
for c = 1:length(colorIDs)
    colors(c,:) = rgb(colorIDs{c});
end

compWeight = {};
allComps = {};

for K = 1:length(models)
    [rxnComps,rxnCompNames,comps{K}] = reactionCompartments(modelStructs{K});
    rxnCompsAnalysis{K} = rxnComps;
    compNames{K} = unique(rxnCompNames);
    allComps = [allComps;comps{K}];
end

allComps = unique(allComps);
colors = colorcube(length(allComps));

modelFields = fieldnames(FBA);

T = 0;
for K = 1:length(models)
    for t = 1:length(Zscore{K}(1,:))+1
        out = FBA.(modelFields{K});
        fluxes = out{t}.x;
        index = [];
        compIndex = [];
        
        for c = 1:length(comps{K})
            weight = sum(abs(fluxes(find(rxnCompsAnalysis{K} == c))));
            compWeight{K}(c,t) = weight;
            compIndex(c,1) = find(ismember(allComps,comps{K}{c}));
        end
    end
end

figure
for K = 1:length(models)
    for t = 1:length(Zscore{K}(1,:))+1
        if t == 1
            compWeight0 = max(compWeight{K},[],2);
        end
        
        compWeight{K}(:,t) = compWeight{K}(:,t)./compWeight0;
        compWeight{K}(isnan(compWeight{K})) = 0;
        
        [~,index] = sort(compWeight{K}(:,t));
        
        
%         labels = [compNames{K}(index)];
        labels = {};
        for i = 1:length(comps{K})
            labels = [labels;{''}];
        end
        
        T = T+1;
        chart = subplot(length(models),length(Zscore{K}(1,:))+1,T);
        pie(compWeight{K}(index,t),labels);
        colormap(chart,colors(compIndex(index),:))
    end
end

%% Compartment Load (molar flux)
figure('visible','on','rend','painters','pos',[10 10 1100 1000])
for K = 1:length(models)
    A0 = max(compWeight{K},[],2);
    
%     allMetCompsTemp = allMetComps{K};
%     ATemp = A{K};cell metabolomics
%     
    A0 = A0(find(~all(compWeight{K}==0,2)));
    allMetCompsTemp = comps{K}(find(~all(compWeight{K}==0,2)));
    ATemp = compWeight{K}(find(~all(compWeight{K}==0,2)),:);
    
    subplot(4,2,K)
    hold on
    box on
    grid on
    for n = 1:length(ATemp(:,1))
        plot(ATemp(n,:)/A0(n),'LineWidth',1.1);
    end
    legend(allMetCompsTemp)
end

%% Compartment element load
colorIDs = {'Black';'LightBlue';'Silver';'Green';'DarkSlateGray';'Pink';'Tomato';'Brown';'PaleGreen';'SaddleBrown';'GoldenRod';'White'};
for c = 1:length(colorIDs)
    colors(c,:) = rgb(colorIDs{c});
end

s = 0;
A = cell(1,length(models));
loadElement = 1; % N4, C1

metaboliteImportance = {};
for K = 1:length(models)
    model = modelStructs{K};
    [~,metComps,comps] = reactionCompartments(model);
    
    [~,~,~,Num] = calculateFormula(model,1:length(model.mets),1);
    
    metaboliteImportance{K} = zeros(length(model.mets),length(Zscore{K}(1,:))+1);
    A{K} = zeros(length(allComps),length(Zscore{K}(1,:))+1);
    
    for r = 1:length(model.rxns)
        compIDInRxn = [];
        rxn = model.rxns{r};
        metsInRxn = findMetsFromRxns(model,rxn);
        metIDsInRxn = findMetIDs(model,metsInRxn);
        compsInRxn = metComps(metIDsInRxn);
        for i = 1:length(compsInRxn)
            compIDInRxn(i,1) = find(ismember(allComps,compsInRxn{i}));
        end
        stoichInRxn = full(model.S(metIDsInRxn,r));
        if length(unique(compsInRxn)) > 1 && ~model.c(r)
            for t = 1:length(Zscore{K}(1,:))+1
                massFlowInRxn = FBA.(models{K}){t}.x(r)*stoichInRxn.*Num(metIDsInRxn,loadElement);
                for m = 1:length(metsInRxn)
                    if massFlowInRxn(m) > 0
                        A{K}(compIDInRxn(m),t) = A{K}(compIDInRxn(m),t) + massFlowInRxn(m);
                        metaboliteImportance{K}(metIDsInRxn(m),t) = metaboliteImportance{K}(metIDsInRxn(m),t) + massFlowInRxn(m);
                    end
                end
            end
        end
    end
end

for K = 1:length(models)
    model = modelStructs{K};
    M = metaboliteImportance{K}(:,6) - metaboliteImportance{K}(:,1);
    M2 = abs(M);
    metImpComp{K} = {};
    [sortM2,indexM] = sort(M2);
    indexM = flip(indexM);
    MinOrder = M(indexM);
    metsInOrder = model.mets(indexM);
    for c = 1:length(allComps)
        i = 1;
        j = 1;
        while i < 100 && j < length(model.mets)
            met = metsInOrder{j};
            weight = MinOrder(j);
            compID = strmatch(met(end),allComps);
            
            if compID == c
                metImpComp{K}{compID}{i,1} = met;
                metImpComp{K}{compID}{i,2} = weight;
                i = i+1;
            end
            j = j+1;
        end
    end
end

figure('visible','on','rend','painters','pos',[10 10 600 800])
for K = 1:length(models)
    A0 = max(A{K},[],2);

    colorPos = find(~all(A{K}==0,2));
    A0 = A0(find(~all(A{K}==0,2)));
    allMetCompsTemp = allComps(find(~all(A{K}==0,2)));
   
    
    ATemp = A{K}(find(~all(A{K}==0,2)),:);
    
    subplot(2,4,K)
    hold on
    box on
    grid on
    for n = 1:length(ATemp(:,1))
        plot(ATemp(n,:)/A0(n),'LineWidth',1.1,'color',colors(colorPos(n),:));
    end
    xlim([1 6])
    ylim([min(min(ATemp./A0))*0.99 1.01])
    legend(allComps(colorPos))
%     set(gca,'YTick',[0.9 0.95 1])
end

%% Compartment activity (mass flow transport reactions)

s = 0;
A = cell(1,length(models));
for K = 1:length(models)
    model = modelStructs{K};
    
    [~,metComps,comps] = reactionCompartments(model);
    
    [~,MW] = calculateFormula(model,1:length(model.mets),1);
    A{K} = zeros(length(allComps),length(Zscore{K}(1,:))+1);
    
    for r = 1:length(model.rxns)
        compIDInRxn = [];
        rxn = model.rxns{r};
        metsInRxn = findMetsFromRxns(model,rxn);
        metIDsInRxn = findMetIDs(model,metsInRxn);
        compsInRxn = metComps(metIDsInRxn);
        for i = 1:length(compsInRxn)
            compIDInRxn(i,1) = find(ismember(allComps,compsInRxn{i}));
        end
        stoichInRxn = full(model.S(metIDsInRxn,r));
        if length(unique(compsInRxn)) > 1 && ~model.c(r)
            for t = 1:length(Zscore{K}(1,:))+1
                massFlowInRxn = FBA.(models{K}){t}.x(r)*stoichInRxn.*MW(metIDsInRxn);
                for m = 1:length(metsInRxn)
                    if massFlowInRxn(m) < 0
                        A{K}(compIDInRxn(m),t) = A{K}(compIDInRxn(m),t) + massFlowInRxn(m);
                    end
                end
            end
        end
    end
end

figure('visible','on','rend','painters','pos',[10 10 1100 1000])
for K = 1:length(models)
    A0 = min(A{K},[],2);
    
%     allMetCompsTemp = allMetComps{K};
%     ATemp = A{K};
    colorPos = find(~all(A{K}==0,2));
    A0 = A0(find(~all(A{K}==0,2)));
    ATemp = A{K}(find(~all(A{K}==0,2)),:);
    
    subplot(4,2,K)
    hold on
    box on
    grid on
    for n = 1:length(ATemp(:,1))
        plot(ATemp(n,:)/A0(n),'LineWidth',1.1,'color',colors(colorPos(n),:));
    end
    legend(allComps(colorPos))
    ylim([0.5 1.05])
end




%% Generate colorbar
N = length(allComps);
rxnCompsAnalysis{K}(end-N+1:end) = 1:N;
comps{K} = allComps;
index = [];
compIndex = [];
T = T+1;
for c = 1:length(comps{K})
    weight = sum(abs(fluxes(find(rxnCompsAnalysis{K} == c))));
    compWeight{K}(c,t) = weight;
    compIndex(c,1) = find(ismember(allComps,comps{K}{c}));
end
index = 1:length(comps{K});

labels = comps{K};

figure
pie(compWeight{K}(index,t),labels);
colormap(colors(compIndex(index),:))
colorbar('Ticks',[1:length(allComps)],'TickLabels',allComps)

%% Metabolite frequency
allMets = {};
for K = 1:length(models)
    allMets = [allMets;metabolites{K}];
end
[uni,~,ic] = unique(allMets);
counts = hist(ic,unique(ic));

maxFreq = max(counts);
maxMets = uni(counts == maxFreq);


%% Compartment heatmap
% colors = {'#ffc700';'#3600ff';'#CC0066';'#4e9800';'#cc0099';'#d7c80f';'#000000';'#404040';'#7edd02';'#cc9999';'#921e1e'};
colorIDs = {'Black';'LightBlue';'Silver';'Green';'DarkSlateGray';'Pink';'Tomato';'Brown';'PaleGreen';'SaddleBrown';'GoldenRod';'White'};
for c = 1:length(colorIDs)
    colors(c,:) = rgb(colorIDs{c});
end

loadFileID = 'BOFsensitivity_';
models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};

DELTA = [];
for K = 1:length(models)
    load([loadFileID,models{K}])
    reactionIDs{K} = rxnID;
    
    for r = 1:length(reactionIDs{K})
        [~,metComps] = reactionCompartments(model);
        compIDInRxn = [];
        rxn = model.rxns{reactionIDs{K}(r)};
        metsInRxn = findMetsFromRxns(model,rxn);
        metIDsInRxn = findMetIDs(model,metsInRxn);
        compsInRxn = metComps(metIDsInRxn);
        for i = 1:length(compsInRxn)
            compIDInRxn(i,1) = find(ismember(allComps,compsInRxn{i}));
        end
        if length(unique(compsInRxn)) > 1 && ~model.c(r)
            rxnComp{K}(r,1) = length(allComps) + 1;
        elseif length(unique(compsInRxn)) == 1 && ~model.c(r)
            rxnComp{K}(r,1) = compIDInRxn(1);
        else
            rxnComp{K}(r,1) = 1;
        end
    end
end

for K = 1:length(models)
    figure
    modelCompartmentHM.(models{K}) = heatmap([rxnComp{K},rxnComp{K}]);
    colormap(colors)
    saveas(gcf,['compartmentHM_',models{K}],'png')
end


save('compartmentHeatmaps.mat')
colorbar
