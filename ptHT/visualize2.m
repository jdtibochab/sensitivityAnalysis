
%% 3D
clear
colors = rand(17,3);

markers = {'o','d','^','s','>','<','p','h'};
models = {'HT','PA','Pt','PtHT','CHO','Yl','Sc'};
% models = {'HT','CHO','Yl','Sc'};
% models = {'PA','Pt'};

figure('visible','on')
hold on
modelP = zeros(length(models),1);
metP = zeros(length(markers),1);
x = [];
y = [];
z = [];
for m = 1:length(models)
    load(['3D_',models{m},'.mat']);
    fieldNames = fieldnames(metTypes);
    groupConectivity = groupConectivity/max(groupConectivity);
    groupCost = groupCost/max(groupCost);
    groupSens = groupSens/max(groupSens);
    
    metTypes.allMets = [];
%     fieldNames = {'AA'};
%     fieldNames = {'TAG'};
%     fieldNames = {'Pe','Pg','Pail','Pchol'};
%     fieldNames = {'As','Dgts','Dgdg','Sqdg','Mgdg'};
    for f = 1:length(fieldNames)
        if ~isempty(metTypes.(fieldNames{f}))
            mets = metTypes.(fieldNames{f});
            for i = 1:length(mets)
                metNum = find(group == mets(i));
                x = [x;groupConectivity(metNum)];
                y = [y;groupCost(metNum)'];
                z = [z;groupSens(metNum)'];
                if groupCod(metNum) > 0
                    markerSize = exp(groupCod(metNum)/max(groupCod))*3;
                    plot3(groupConectivity(metNum),groupCost(metNum),groupSens(metNum),markers{m},'MarkerFace',colors(f,:),'MarkerEdge',colors(f,:),'MarkerSize',markerSize);
                    
                end
            end
            if m == 1
                metP(f) = plot(0,0,'s','MarkerFace',colors(f,:),'MarkerEdge',colors(f,:));
            end
        end
    end
    modelP(m) = plot(0,0,markers{m},'MarkerFace','k','MarkerEdge','k');
end
x(find(isnan(x))) = 0;
y(find(isnan(y))) = 0;
z(find(isnan(z))) = 0;

poly = 'poly11';

surffit1 = fit([x,y],z,poly);
fdata = feval(surffit1,[x,y]);
e = abs(fdata - z) > 1.5*std(z);

surffit = fit([x,y],z,poly,'Exclude',e);
confint(surffit)';
% outliers = excludedata(y,z,'indices',e);

plot(surffit)

xlim([0 1])
ylim([0 1])
zlim([0 1])

legend(modelP,models)

xlabel('Conectivity')
ylabel('Cost')
zlabel('Sensitivity')

%%
saveas(gcf,'3Dplot_ALL_RS','svg')

%% Parameter analysis
clear
colors = rand(17,3);

modelGroups.all = {'HT','PA','Pt','CHO','Yl','Sc'};
modelGroups.ht = {'HT','CHO','Yl','Sc'};
modelGroups.pa = {'PA','Pt'};

metGroups.all = {'allMets'};
metGroups.AA = {'AA'};
metGroups.TAG = {'TAG'};
metGroups.P = {'Pe','Pg','Pail','Pchol'};
metGroups.L = {'As','Dgts','Dgdg','Sqdg','Mgdg'};

modelFieldNames = fieldnames(modelGroups);
metFieldNames = fieldnames(metGroups);
P = [];
for g = 1:length(modelFieldNames)
    models = modelGroups.(modelFieldNames{g});
    Gi(g) = g;
    for n = 1:length(metFieldNames)
        MGi(g) = n;
        fieldNames = metGroups.(metFieldNames{n});
        N = 0;
        x = [];
        y = [];
        z = [];
        for m = 1:length(models)
            load(['3D_',models{m},'.mat']);
            groupConectivity = groupConectivity;%/max(groupConectivity);
            groupCost = groupCost;%/max(groupCost);
            groupSens = groupSens;%/max(groupSens);
            for f = 1:length(fieldNames)
                if ~isempty(metTypes.(fieldNames{f}))
                    mets = metTypes.(fieldNames{f});
                    for i = 1:length(mets)
                        metNum = find(group == mets(i));
                        x = [x;groupConectivity(metNum)];
                        y = [y;groupCost(metNum)'];
                        z = [z;groupSens(metNum)'];
                        N = N + length(metNum);
                    end
                end
            end
        end
        x(find(isnan(x))) = 0;
        y(find(isnan(y))) = 0;
        z(find(isnan(z))) = 0;
        
%         [x,index] = sort(x);
%         z = z(index);
        
%         poly = 'poly11';
%         surffit1 = fit([x,y],z,poly);
%         fdata = feval(surffit1,[x,y]);
%         e = abs(fdata - z) > 1.5*std(z);
%         [surffit,gof] = fit([x,y],z,poly);

        poly = 'poly2';
        [surffit,gof] = fit(x,z,poly);
        fdata = feval(surffit,x);
        figure;
        plot(x,fdata,x,z,'.');
        
        parVal = coeffvalues(surffit)';
        confInt = confint(surffit)';
        S = (confInt(:,2)-confInt(:,1))/3.92;
        R2 = gof.adjrsquare;
        if R2<0
            R2 = 0;
        end
        t = sqrt(R2)*sqrt((N-2)/(1-R2));
        corrPval = 1-tcdf(t,N-2);
        
        gi =  ones(length(parVal),1)*g;
        mgi = ones(length(parVal),1)*n;
        pi = [1:length(parVal)]';
        Ni = ones(length(parVal),1)*N;
        R2i = ones(length(parVal),1)*R2;
        ti = ones(length(parVal),1)*t;
        corrPvali = ones(length(parVal),1)*corrPval;
        
        Pi = [gi,mgi,pi,parVal,S,Ni,R2i,ti,corrPvali];
        
        P = [P;Pi];
    end
end

% metGroup comparison
for g = 1:length(modelFieldNames)
    rowsG = find(P(:,1) == Gi(g));
    parMat = P(rowsG,:);
    for p = 1:length(pi)
        rowsP = find(parMat(:,3) == pi(p));
        parMati = parMat(rowsP,:);
        for i = 1:length(rowsP)
            for j = 1:length(rowsP)
                if i~=j
                    M1 = parMati(i,4);
                    M2 = parMati(j,4);
                    S1 = parMati(i,5);
                    S2 = parMati(j,5);
                    N1 = parMati(i,6);
                    N2 = parMati(j,6);
                    [p_temp,t_temp] = statDiff(M1,M2,S1,S2,N1,N2);
                    pval.(modelFieldNames{g}).(['par',num2str(pi(p))])(i,j) = p_temp;
                    tval.(modelFieldNames{g}).(['par',num2str(pi(p))])(i,j) = t_temp;
                end
            end
        end
    end
end

%% Parameter analysis with means
clear
colors = rand(17,3);

modelGroups.all = {'HT','PA' ,'Pt','CHO','Yl','Sc'};
modelGroups.ht = {'HT','CHO','Yl','Sc'};
modelGroups.pa = {'PA','Pt'};

metGroups.all = {'AA','TAG','Pe','Pg','Pail','Pchol','As','Dgts','Dgdg','Sqdg','Mgdg'};
% metGroups.AA = {'AA'};
% metGroups.TAG = {'TAG'};
% metGroups.P = {'Pe','Pg','Pail','Pchol'};
% metGroups.L = {'As','Dgts','Dgdg','Sqdg','Mgdg'};

modelFieldNames = fieldnames(modelGroups);
metFieldNames = fieldnames(metGroups);
P = [];
SP = 0;
figure;
for g = 1:length(modelFieldNames)
    models = modelGroups.(modelFieldNames{g});
    Gi(g) = g;
    for n = 1:length(metFieldNames)
        MGi(g) = n;
        fieldNames = metGroups.(metFieldNames{n});
        N = 0;
        x = [];
        y = [];
        z = [];
        gCon = 0;
        gCost = 0;
        gSens = 0;
        for m = 1:length(models)
            load(['3D_',models{m},'.mat']);
            groupConectivity = groupConectivity/max(groupConectivity);
            groupCost = groupCost/max(groupCost);
            groupSens = groupSens/max(groupSens);
            
            for f = 1:length(fieldNames)

                if ~isempty(metTypes.(fieldNames{f}))
                    mets = metTypes.(fieldNames{f});
                    metNum = find(ismember(group,mets));
                    
                    gCon(m,f) = mean(groupConectivity(metNum));
                    gCost(m,f) = mean(groupCost(metNum));
                    gSens(m,f) = mean(groupSens(metNum));
                end
            end
        end
        
        for i = 1:length(fieldNames)            
            temp = gCon(:,i);
            x(i,1) = mean(temp(temp~=0));
            stdX(i,1) = std(temp(temp~=0));

            temp = gCost(:,i);
            y(i,1) = mean(temp(temp~=0));
            stdY(i,1) = std(temp(temp~=0));
            
            temp = gSens(:,i);
            z(i,1) = mean(temp(temp~=0));
            stdZ(i,1) = std(temp(temp~=0));
        end
        x(isnan(z)) = [];
        stdX(isnan(z)) = [];
        y(isnan(z)) = [];
        stdY(isnan(z)) = [];
        stdZ(isnan(z)) = [];
        fieldNames(isnan(z)) = [];
        
        z(isnan(z)) = [];
        [x,index] = sort(x);
        y = y(index);
        z = z(index);
        fieldNames = fieldNames(index);
        
        
        poly = 'poly1';
        surffit = fit(y,z,poly);
        fdataY = feval(surffit,y);
        surffit2 = fit(x,z,poly);
        fdataX = feval(surffit2,x);
        
        SP = SP + 1;
        subplot(2,length(modelFieldNames),SP);
        hold on
        plot(y,fdataY);
        errorbar(y,z,stdZ,stdZ,stdY,stdY,'.');
        text(y,z,fieldNames)
        xlabel('Cost')
        ylabel('Sensitivity')
        
        subplot(2,length(modelFieldNames),SP+3);
        hold on
        plot(x,fdataX);
        errorbar(x,z,stdZ,stdZ,stdX,stdX,'.');
        text(x,z,fieldNames)
        xlabel('Conectivity')
        ylabel('Sensitivity')
        
    end
end
%%
saveas(gcf,'6plotCostConec','svg')
%% Box plots
clear
colors = rand(15,3);

markers = {'o','d','^','s','>','<','p','h'};
models = {'HT','PA','Pt','CHO','Yl','Sc'};%,'PtHT'};

models = {'CHO','Yl','Sc','HT'};
% models = {'Pt','PA'};

groupSensBP = {};

i = 0;
for m = 1:length(models)

    load(['3D_',models{m},'.mat']);
    fieldNames = fieldnames(metTypes);
    metTypes.allMets = [];
    metTypes.rest = [];
    for f = 1:length(fieldNames)
        if f > length(groupSensBP)
            groupSensBP{f} = [];
        end
        if ~isempty(metTypes.(fieldNames{f}))
            mets = metTypes.(fieldNames{f});
            metNum = find(ismember(group,mets));
            data = growth_sens_vector(metNum,:)/max(growth_sens_vector(:));
            groupSensBP{f} = [groupSensBP{f};data(:)];
        end
    end
end

for f = 1:length(fieldNames)
    medians(f) = median(groupSensBP{f});
end

[~,index] = sort(medians);

groupSensVec =[];
g = [];
for i = 1:length(fieldNames)
    f = index(i);
    groupSensVec = [groupSensVec;groupSensBP{f}];
    g = [g;(i-1)*ones(length(groupSensBP{f}),1)];
end

figure('visible','on');
boxplot(groupSensVec,g);
set(gca, 'XTickLabel', fieldNames(index))
ylabel('Sensitivity')
% title('Autotrophy')
title('Heterotrophy');
%%
saveas(gcf,'HetBoxPlotAA','svg')
% saveas(gcf,'PABoxPlotAA','svg')
%% Sample matrix
clear

models = {'HT','PA','Pt','CHO','Yl','Sc'};
modelsPA = {'PA','Pt',[],[],[]};
modelsHT = {'HT','PtHT','CHO','Yl','Sc'};

models = [modelsHT',modelsPA'];

dim = size(models);
N1 = dim(1);
N2 = dim(2);

height = 0.8 / N1;%%
saveas(gcf,'6plotCostConec','svg')
width = 0.8 / N2;

figure('visible','on')

for m = 1:N1
    for n = 1:N2
        if ~isempty(models{m,n})
            left = 0.1 + (n - 1) * width;
            bottom = 0.9 - m * height;
            
            load(['3D_',models{m,n},'.mat']);
            fieldNames = fieldnames(metTypes);
            
            groupConectivity = groupConectivity/max(groupConectivity);
            groupCost = groupCost/max(groupCost);
            groupSens = groupSens/max(groupSens);
            groupCod = groupCod/max(groupCod);
            
%             fieldNames = {'AA'};
            subplot('position', [left bottom width height]);
            hold on
            for f = 1:length(fieldNames)
                if ~isempty(metTypes.(fieldNames{f}))
                    mets = metTypes.(fieldNames{f});
                    for i = 1:length(mets)
                        metNum = find(group == mets(i));
                        plot(groupCost(metNum),groupSens(metNum),'.k');
                    end
                end
            end
            ylim([0 1])
            xlim([0 1])
            set(gca, 'YTickLabel', []);
            set(gca, 'XTickLabel', []);
            set(gca, 'YTick', []);
            set(gca, 'XTick', []);
        end
    end
end




%% ALL (AA)

clear
models = {'HT','PA','Pt','Yl','Sc','CHO'};

S = [];
for m = 1:length(models)
    load(['3D_',models{m},'.mat']);
    Mets{m} = Met(group);
    modelSens{m} = groupSens/max(groupSens);
    modelCod{m} = groupCod/max(groupCod);
end

load(['3D_PA.mat']);
intMet = Met(metTypes.AA);

for m = 1:length(models)
    metNum{m} = find(ismember(intMet,Mets{m}));
end

groupSens = zeros(length(metNum{1}),length(models));

for m = 1:length(models)
    met = intMet(metNum{m});
    num = find(ismember(Mets{m},met));
    for e = 1:length(intMet)
        index = strmatch(intMet{e},Mets{m},'exact');
        if ~isempty(index)
            S(e,m) = modelSens{m}(index);
            C(e,m) = modelCod{m}(index);
        else
            S(e,m) = 0;
            C(e,m) = 0;
        end
    end
end

C(:,end+1) = sum(C(:,1:end),2);
[C,indexC] = sortrows(C,length(C(1,:)));

S(:,end+1) = sum(S(:,1:end),2);
[S,indexS] = sortrows(S,length(S(1,:)));

textMet = intMet;

for i = 1:length(textMet)
    met = char(textMet{i});
    a = met(1);
    a = upper(a);
    met(1) = a;
    
    met(end-1:end) = [];
    eliminate = strfind(met,'__L');
    met(eliminate:eliminate+2) = [];
    
    textMet{i} = met;
end

figure
bar(S(:,1:end-1),'stacked')
set(gca, 'XTick', [1:length(intMet)])
set(gca, 'XTickLabel', textMet(indexS))
xlabel('Metabolites')
ylabel('Sensitivity')
legend(models,'Location','NorthWest')
saveas(gcf,'AAmets_S_ALL','svg')

figure
bar(C(:,1:end-1),'stacked')
set(gca, 'XTick', [1:length(intMet)])
set(gca, 'XTickLabel', textMet(indexC))
xlabel('Metabolites')
ylabel('Codependence')
legend(models,'Location','NorthWest')

saveas(gcf,'AAmets_C_ALL','svg')

%% ALL
clear
models = {'HT','PA','Pt','CHO','Yl','Sc'};

S = [];
for m = 1:length(models)
    load(['3D_',models{m},'.mat']);
    Mets{m} = Met(group);
    modelSens{m} = groupSens;
    modelCod{m} = groupCod;
end

allMetsFromModels = {};
for m = 1:length(models)
    allMetsFromModels = [allMetsFromModels;Mets{m}];
end
allMetsFromModels = unique(allMetsFromModels);

for e = 1:length(allMetsFromModels)
    for m = 1:length(models)
        metNum = find(ismember(allMetsFromModels{e},Mets{m}));
        if ~isempty(metNum)
            S(e,m) = modelSens{m}(metNum);
            C(e,m) = modelCod{m}(metNum);
        else
            S(e,m) = 1;
            C(e,m) = 1;
        end
    end
end
C = -log2(C);
C(:,end+1) = sum(C(:,1:end),2);
[C,indexC] = sortrows(C,length(C(1,:)));

S = -log2(S);
S(:,end+1) = sum(S(:,1:end),2);
[S,indexS] = sortrows(S,length(S(1,:)));

figure
bar(S(:,1:end-1),'stacked')
xlabel('Metabolites')
ylabel('Sensitivity')
legend(models,'Location','NorthWest')
saveas(gcf,'Barplot_S_ALL','svg')

figure
bar(C(:,1:end-1),'stacked')
xlabel('Metabolites')
ylabel('Codependence')
legend(models,'Location','NorthWest')

saveas(gcf,'Barplot_C_ALL','svg')

%% All, lines
clear
models = {'HT','PA','Pt','CHO','Yl','Sc'};
% models = {'Sc'}
S = [];
for m = 1:length(models)
    load(['3D_',models{m},'.mat']);
    Mets{m} = Met(group);
    modelSens{m} = groupSens;%/max(Sens);
    modelCod{m} = groupCod;
end

allMetsFromModels = {};
for m = 1:length(models)
    allMetsFromModels = [allMetsFromModels;Mets{m}];
end
allMetsFromModels = unique(allMetsFromModels);

figure
subplot(1,2,1)
hold on
for m = 1:length(models)
    plot(1:numel(modelSens{m}),sort(modelSens{m}),'o')
end
xlabel('Metabolites')
ylabel('Sensitivity')
legend(models)

subplot(1,2,2)
hold on
for m = 1:length(models)
    plot(1:numel(modelCod{m}),sort(modelCod{m}),'o')
end
legend(models)
xlabel('Metabolites')
ylabel('Codependence')
saveas(gcf,'lines_C_ALL','svg')
