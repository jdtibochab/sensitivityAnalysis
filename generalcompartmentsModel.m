function [G_subs,G_subs_corrected] = generalcompartmentsModel(subSystems)
%Cristal Zuniga
%a=iSA_NCTC.subSystems;

a=subSystems;

% a=model.subSystems;
Subs=unique(a);
%Subs=a;

N_subs={};
for i=1:numel(Subs)
    index = find(ismember(Subs{i}, ';'));
    str=Subs{i};
    if ~isempty(index)
        if length(index)==1
            [token, remain] = strtok(str, ';');
            N_subs{end+1,1}=token;
            N_subs{end+1,1}=remain(2:end);
        else
            [token, remain] = strtok(str, ';');
            N_subs{end+1,1}=token;
            [token, remain]= strtok(remain(2:end), ';');
            N_subs{end+1,1}=token;
            N_subs{end+1,1}=remain(2:end);
        end
    else
        N_subs{end+1,1}=str;
    end
end

N_subs2=unique(N_subs);

for i=1:numel(N_subs2)
    IndexC = strfind(N_subs,N_subs2{i});
    Index = find(not(cellfun('isempty', IndexC)));
    N_subs2{i,2}=length(Index);
end


%Looking for unique names and counting

for i=1:length(N_subs2)
    A=strmatch(N_subs2{i},a,'exact');  %
    N_subs2{i,3}=length(A);
end

%Check should be of the same size of model.rxns
Check=sum(cellfun(@double,N_subs2(:,3)));

%Group by category
AAs={'alanine','aspartate','arginine','proline','cystein','methionine','quinate',...
    'glutamate','glycine','threonine','histidine','lysine','phenylalanine','methylglyoxal',...
    'tyrosine','tryptophan','valine','leucine','isoleucine','asparagine','glutamine','serine','trna','trna charging','protein',...
    'amino acidand','hydroxyphenylacetate','translation','amino acid','rna','dehydrogenases','hydrogenases','ribosome','auxin','chorismate'};
Faty={'triacylglycerol','r group','lipid','phosphatidylethanolamine','fatty','glycerolipid','glycerophospholipid','sphingolipid','inositol','sulfolipid','sterol',...
    'linoleic','phospholipid','alpha-linolenic','arachidonic','cardiolipin',...
    'teichoic'};
CofVit={'vitamin','thiamine','riboflavine','riboflavin','vitamin b6','nicotinate','pantothenate','biotin','lipoic','folate','thioredoxin',...
    'porphyrin','chlorophyll','ubiquinone','retinol','ascorbate','cofactor','pyridoxine','nad','selenoamino','coenzymeand',...
    'thiamin','coenzyme','siroheme','cobalamin'};
Glycan={};
OAAs={'beta-alanine','taurine','selecompound','cyanoamino','glutathione','cyanophycin','aminosugars','other','phosphinate','aminosugar'};
terpenoid={'bile','terpenoid','carotenoid','carotene','isoprenoid','quinone','brassinosteroid','polyamine','steroids','steroid','polyprenyl','putrescine'};
Xen={'benzoate','aminobenzoate','fluorobenzoate','toluene','xylene','nitrotoluene','ddt','bispherol','dioxin'...
    'naphthalene','polycyclic','drug','xenobiotic','aromatic'};
Carb={'glycan','glycosaminoglycan','glycosylphosphatidylinositol','galactolipids','glycosphingolipid','lipopolysaccharide','glycoprotein','murein','glycolysis','gluconeogenesis','c5-branched','acetoin','butanediol',...
    'citrate cycle','tca cycle','citric acid cycle','anaplerotic'...
    'pentose phosphate pathway',...
    'glucuronate','galactarate','glucarate','gluconate',...
    'fructose','mannose','rhamnose',...
    'alactose','glycan','sialic',...
    'aldarate','glycerate','muconate',...
    'starch', 'sucrose','carbon','fermentations',...
    'amino sugar and nucleotide sugar',...
    'pyruvate','glyoxylate','dicarboxylate',...
    'propanoate','ketone bodies','exopolysaccharide',...
    'butanoate', 'arabinose','alternate carbon', 'xylose','alcohol','sugar','carbohydrateand','alkane','alkanesulfonate','cmp'};
Energy={'oxidative_phosphorylation','citric_acid','oxidative phosphorylation',...
    'photosynthesis','chemotaxis',...
    'antenna','nitrosative',...
    'fixation',...
    'methane',...
    'nitrogen',...
    'sulfur','urea','atpase','atp','ammonia','nitrate'};
Nucleotide={'nucleotide','purine','pyrimidine','nucleosidases','deoxyribose','deoxynucleoside','nucleoside','ppgpp','xanthine'};
Biomass={'biomass','maintenance','demand','spectral','light','diphtamide','osmoregulation','tungstate','pilus','phosphate','potassium'};
CellE={'cell','envelope','wall','membrane'};
Transport = {'transport'};

G_subs_corrected = a;
G_subs=a;Aminoacids=0;Fattys=0;Transport=0;cofVit=0;a=0;b=0;c=0;d=0;e=0;f=0;gly=0;oaa=0;ter=0;xeno=0;car=0;g=0;h=0;m=0;ene=0;nuc=0;n=0;bio=0;x=0;cell=0;
for i=1:length(G_subs)
    a=0;b=0;c=0;d=0;e=0;f=0;g=0;h=0;m=0;n=0;x=0;cell=0;
    found = 0;
    A=lower(G_subs{i});
    for j=1:numel(OAAs)
        k = strfind(A, OAAs{j});
        if ~isempty(k)
            %oaa=oaa+G_subs{i,3};
            G_subs{i,1}='Metabolism of other amino acids';
            a=1;
            found = 1;
            break
        end
    end
    if a~=1
        for j=1:numel(AAs)
            k = strfind(A, AAs{j});
            if ~isempty(k)
                %Aminoacids=Aminoacids+G_subs{i,3};
                G_subs{i,1}='Amino acid metabolism';
                b=1;
                found = 1;
                break
            end
        end
        
        if b~=1
            for j=1:numel(Faty)
                k = strfind(A, Faty{j});
                if ~isempty(k)
                    %Fattys=Fattys+G_subs{i,3};
                    G_subs{i,1}='Fatty acid metabolism';
                    c=1;
                    found = 1;
                    break
                end
            end
            if c~=1
                for j=1:numel(CofVit)
                    k = strfind(A, CofVit{j});
                    if ~isempty(k)
                        %cofVit=cofVit+G_subs{i,3};
                        G_subs{i,1}='Metabolism of cofactors and vitamins';
                        d=1;
                        found = 1;
                        break
                    end
                end
                if d~=1
                    for j=1:numel(Glycan)
                        k = strfind(A, Glycan{j});
                        if ~isempty(k)
                            %gly=gly+G_subs{i,3};
                            G_subs{i,1}='Glycan biosynthesis and metabolism';
                            e=1;
                            found = 1;
                            break
                        end
                    end
                    if e~=1
                        for j=1:numel(terpenoid)
                            k = strfind(A, terpenoid{j});
                            if ~isempty(k)
                                %ter=ter+G_subs{i,3};
                                G_subs{i,1}='Terpenoid biosynthesis';
                                f=1;
                                found = 1;
                                break
                            end
                        end
                        if f~=1
                            for j=1:numel(Xen)
                                k = strfind(A, Xen{j});
                                if ~isempty(k)
                                    % xeno=xeno+G_subs{i,3};
                                    G_subs{i,1}='Xenobiotics biodegradation and metabolism';
                                    g=1;
                                    found = 1;
                                    break
                                end
                            end
                            if g~=1
                                for j=1:numel(Carb)
                                    k = strfind(A, Carb{j});
                                    if ~isempty(k)
                                        %car=car+G_subs{i,3};
                                        G_subs{i,1}='Carbohydrate metabolism';
                                        h=1;
                                        found = 1;
                                        break
                                    end
                                end
                                if h~=1
                                    for j=1:numel(Energy)
                                        k = strfind(A, Energy{j});
                                        if ~isempty(k)
                                            %ene=ene+G_subs{i,3};
                                            G_subs{i,1}='Energy metabolism';
                                            n=1;
                                            found = 1;
                                            break
                                        end
                                    end
                                    if n~=1
                                        for j=1:numel(Nucleotide)
                                            k = strfind(A, Nucleotide{j});
                                            if ~isempty(k)
                                                %nuc=nuc+G_subs{i,3};
                                                G_subs{i,1}='Nucleotide metabolism';
                                                x=1;
                                                found = 1;
                                                break
                                            end
                                        end
                                        if x~=1
                                            for j=1:numel(Biomass)
                                                k = strfind(A, Biomass{j});
                                                if ~isempty(k)
                                                    %bio=bio+G_subs{i,3};
                                                    G_subs{i,1}='Biomass, maintenance, demand, and expectral decomposition';
                                                    cell=1;
                                                    found = 1;
                                                    break
                                                end
                                            end
                                            if cell~=1
                                                for j=1:numel(CellE)
                                                    k = strfind(A, CellE{j});
                                                    if ~isempty(k)
                                                        %cell=cell+G_subs{i,3};
                                                        G_subs{i,1}='Cell envelope';
                                                        m=1;
                                                        found = 1;
                                                        break
                                                    end
                                                end
                                                
                                                
                                                if m~=1
                                                    k = strfind(A, 'transport');
                                                    kkk = strfind(A, 'exchange');
                                                    if ~isempty(k) || ~isempty(kkk)
                                                        %Transport=Transport+G_subs{i,3};
                                                        G_subs{i,1}='Transport reactions';
                                                        
                                                        found = 1;
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        
    end
    G_subs_corrected{i,1} = G_subs{i,1};
    if ~found
        G_subs_corrected{i,1} = 'Other';
    end
end