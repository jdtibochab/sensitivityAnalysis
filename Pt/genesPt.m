%% Genes in CCM

load('iLB1027_lipid.mat')
model = iLB1027_lipid;
model.subSystems;

central = {'carbohydrate','glycolysis', 'TCA','Calvin','Pyruvate', 'carbon fixation', 'pentose phosphate'};

ss = [];
inc = [];
exc = [];
allGenes = model.genes;
CCMgenes = [];
for r = 1:length(model.rxns)
    c = 0;
    for i = 1:length(central)
        if contains(model.subSystems{r},central(i),'IgnoreCase',true)
            ss = [ss;model.subSystems(r)];  
            c = 1;
        end
    end
    if ~c
        exc = [exc;model.subSystems(r)];
    else
        inc = [inc;model.subSystems(r)]; 
        reactionGeneIDs = findGenesFromRxns(model,model.rxns{r});
        reactionGeneIDs{1}
        CCMgenes = [CCMgenes;reactionGeneIDs{1}];
    end
end
exc = unique(exc)
inc = unique(inc)
CCMgenes = unique(CCMgenes);

%%
geneIDs = findGeneIDs(model,CCMgenes);
for i = 1 :length(CCMgenes)
    rxnIDs = find(model.rxnGeneMat(:,geneIDs(i)));
    model.subSystems{rxnIDs}
end

%%
%% Genes in CCM

load('iLB1027_lipid.mat')
model = iLB1027_lipid;
model.subSystems;

central = {'lipid'};

ss = [];
inc = [];
exc = [];
allGenes = model.genes;
CCMgenes = [];
for r = 1:length(model.rxns)
    c = 0;
    for i = 1:length(central)
        if contains(model.subSystems{r},central(i),'IgnoreCase',true)
            printRxnFormula(model,model.rxns{r});
        end
    end
end
