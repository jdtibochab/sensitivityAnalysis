function [metTypes,metGroups,shortNames] = findMetType(model,mets)

Met_ID = findMetIDs(model,mets);
shortNames = mets;
% Lipids
TAG = strmatch('tag',mets);
    TAG = [TAG;strmatch('m15',mets);strmatch('triglyc',mets)];
DAG = strmatch('dag',mets);
MAG = strmatch('mag',mets);
    MAG = [MAG;strmatch('pmtcoa',mets);strmatch('stcoa',mets)];
Pe = strmatch('pe',mets);
Pg = strmatch('pg',mets);
As = strmatch('asqd',mets);
Dgts = strmatch('dgts',mets);
Dgdg = strmatch('dgdg',mets);
Sqdg = strmatch('sqdg',mets);
Mgdg = strmatch('mgdg',mets);
Pail = strmatch('pail',mets);
Pchol = strmatch('pchol',mets);

for i=1:length(TAG)
    shortNames(TAG(i))={['tag',num2str(i)]};
end
for i=1:length(DAG)
    shortNames(DAG(i))={['dag',num2str(i)]};
end
for i=1:length(MAG)
    shortNames(MAG(i))={['mag',num2str(i)]};
end
for i=1:length(Dgts)
    shortNames(Dgts(i))={['dgts',num2str(i)]};
end
for i=1:length(Dgdg)
    shortNames(Dgdg(i))={['dgdg',num2str(i)]};
end
for i=1:length(As)
    shortNames(As(i))={['asqd',num2str(i)]};
end
for i=1:length(Pg)
    shortNames(Pg(i))={['pg',num2str(i)]};
end
for i=1:length(Pe)
    shortNames(Pe(i))={['pe',num2str(i)]};
end
for i=1:length(Sqdg)
    shortNames(Sqdg(i))={['sqdg',num2str(i)]};
end
for i=1:length(Mgdg)
    shortNames(Mgdg(i))={['mgdg',num2str(i)]};
end

    
% AA
AA = [];
AA1 = strfind(mets,'__L');
AA2 = strfind(mets,'gly_');

j = 1;
for i = 1:length(AA1)
   if AA2{i} == 1
       AA(j,1) = i;
       j = j+1;
   end
   if AA1{i} == 4
       AA(j,1) = i;
       j = j+1;
   end
end

AA(strmatch('gal__L_c',mets(AA))) = [];
AA(strmatch('mal__L_c',mets(AA))) = [];

for i = 1:length(AA)
   shortNames(AA(i)) = strtok(mets(AA(i)),'_'); 
   a = char(shortNames(AA(i)));
   b = upper(a(1));
   txt = [b,a(2:end)];
   shortNames(AA(i)) = {txt};
end

% CB
CB = [];
CB1 = strfind(model.metNames(Met_ID),'ose');
notCB1 = strfind(model.metNames(Met_ID),'serine');
j = 1;
for i = 1:length(CB1) 
   if CB1{i} ~= 0
       if notCB1{i} ~= 0
       else
           CB(j,1) = i;
           j = j+1;
       end
   end
end
CB = [CB;strmatch('lac__D',mets)];

% Nuc
Nuc = [];
Nuc1 = strfind(mets,'tp');
j = 1;
for i = 1:length(Nuc1)
   if Nuc1{i} ~= 0
       Nuc(j,1) = i;
       j = j+1;
   end
end

metTypes.TAG = unique(TAG);
metTypes.DAG = unique(DAG);
metTypes.MAG = unique(MAG);
metTypes.Pe = unique(Pe);
metTypes.Pg = unique(Pg);
metTypes.As = unique(As);
metTypes.Dgts = unique(Dgts);
metTypes.Dgdg = unique(Dgdg);
metTypes.Sqdg = unique(Sqdg);
metTypes.Mgdg = unique(Mgdg);
metTypes.Pail = unique(Pail);
metTypes.Pchol = unique(Pchol);

metTypes.AA = unique(AA);

metTypes.CB = unique(CB);

metTypes.Nuc = unique(Nuc);

fieldNames = fieldnames(metTypes);
allMets = [];
for  f = 1:length(fieldNames)
    allMets = [allMets;metTypes.(fieldNames{f})];
end

all = 1:length(mets);
metTypes.rest = setdiff(all',allMets);

allMets = [allMets;metTypes.rest];
allMets = unique(allMets);

metTypes.allMets = allMets;

%%
metGroups.AA = {'AA'};
metGroups.CB = {'CB'};
metGroups.Nuc = {'Nuc'};
metGroups.FA = {'TAG','DAG','MAG'};
metGroups.PL = {'Pe','Pg','Pail','Pchol'};
metGroups.OL = {'As','Dgts','Dgdg','Sqdg','Mgdg'};
metGroups.rest = {'rest'};
end