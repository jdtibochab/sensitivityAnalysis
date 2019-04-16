%% Generate BOF table from .csv files
load('iCHOv1.mat')
model = iCHOv1;

myDir = dir('.');
filenames={myDir(:).name}';
files = filenames(endsWith(filenames,'.csv'));
mets = strtok(files,'.');

for m=1:length(mets)
    mets{m} = [lower(mets{m}),'_c'];
    metIDs(m) = findMetIDs(model,mets{m});
    if metIDs(m)
        foundMets{m,1} = model.mets(metIDs(m));
    end
end

opts = detectImportOptions(files{1});

for f = 1:length(files)
    T = readtable(files{f},opts);
    Stoich(f,:) = T.Var2; 
end
time = (round(T.Var1))';
M = [time;Stoich];