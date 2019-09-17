
%% 3D
clear
colors = rand(17,3);

markers = {'o','d','^','s','>','<','p','h'};
models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};
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
            metNums = metTypes.(fieldNames{f});
            for i = 1:length(metNums)
                metNum = find(group == metNums(i));
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

xlabel('Connectivity')
ylabel('Cost')
zlabel('Sensitivity')

%%
saveas(gcf,'3Dplot_ALL_RS','svg')

%% Parameter analysis
clear
colorIDs = {'#ff7f00';'#e41a1c';'#984ea3';'#ffffb3';'#377eb8';'#4daf4a';'#a65628'};
colors = hex2rgb(colorIDs);

modelGroups.all = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};
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
                    metNums = metTypes.(fieldNames{f});
                    for i = 1:length(metNums)
                        metNum = find(group == metNums(i));
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
colorIDs = {'#ff7f00';'#e41a1c';'#984ea3';'#ffffb3';'#377eb8';'#4daf4a';'#a65628'};
colors = hex2rgb(colorIDs);


models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};
modelGroups.ht = {'HT','CHO','Yl','Sc','PtHT'};
modelGroups.pa = {'PA','Pt'};

metGroups.all = {'CB','AA','TAG','Pe','Pg','Pail','Pchol','As','Dgts','Dgdg','Sqdg','Mgdg'};
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
        gCod = 0;
        stdX = 0;
        stdY = 0;
        for m = 1:length(models)
            load(['3D_',models{m},'.mat']);
            groupConectivity = groupConectivity/max(groupConectivity);
            groupCost = groupCost/max(groupCost);
            groupSens = groupSens/max(groupSens);
            
            for f = 1:length(fieldNames)

                if ~isempty(metTypes.(fieldNames{f}))
                    metNums = metTypes.(fieldNames{f});
                    metNum = find(ismember(group,metNums));
                    
                    gCon(m,f) = mean(groupConectivity(metNum));
                    gCod(m,f) = mean(groupCod(metNum));
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
            
%             temp = gSens(:,i);
            temp = gCod(:,i);
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
        [surffit,gof] = fit(y,z,poly);
        fdataY = feval(surffit,y);
        [surffit2,gof2] = fit(x,z,poly);
        fdataX = feval(surffit2,x);
        
        gof.rsquare
        gof2.rsquare
        
        
        SP = SP + 1;
        subplot(2,length(modelFieldNames),SP);
        hold on
        grid on
        box on
        plot(y,fdataY);
        errorbar(y,z,stdZ,stdZ,stdY,stdY,'.');
        text(y,z,fieldNames)
        xlabel('Relative Cost')
%         ylabel('Relative Sensitivity')
        ylabel('Codependence')
        
        colormap(colorcube)
        
        subplot(2,length(modelFieldNames),SP+3);
        hold on
        grid on
        box on
        plot(x,fdataX);
        errorbar(x,z,stdZ,stdZ,stdX,stdX,'.');
        text(x,z,fieldNames)
        xlabel('Relative Connectivity')
%         ylabel('Relative Sensitivity')
        ylabel('Codependence')
        
        colormap(colorcube)
        
        
    end
end
%%
% saveas(gcf,'6plotSensitivity','svg')
saveas(gcf,'6plotCodependence','svg')
%% Box plots
clear
colors = rand(15,3);

markers = {'o','d','^','s','>','<','p','h'};
% models = {'HT','PA','Pt','PtHT','CHO','Yl','Sc'};

% models = {'CHO','Yl','Sc','HT','PtHT'};
models = {'Pt','PA'};

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
            metNums = metTypes.(fieldNames{f});
            metNum = find(ismember(group,metNums));
            data = groupSens(metNum)/max(groupSens);
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
title('Autotrophy')
% title('Heterotrophy');
%%
% saveas(gcf,'HetBoxPlotAA','svg')
saveas(gcf,'PABoxPlotAA','svg')

%% Sample matrix
clear

models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};
modelsPA = {'PA','Pt',[],[],[]};
modelsHT = {'HT','PtHT','CHO','Yl','Sc'};

models = [modelsHT',modelsPA'];

dim = size(models);
N1 = dim(1);
N2 = dim(2);

height = 0.8 / N1;
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
                    metNums = metTypes.(fieldNames{f});
                    for i = 1:length(metNums)
                        metNum = find(group == metNums(i));
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

%%
saveas(gcf,'sampleMatrix','svg')

%% ALL (AA)

clear




models = {'PA','HT','Pt','PtHT','Sc','Yl','CHO'}; G = 'ALL';colorIDs = {'#4daf4a';'#a6d854';'#8da0cb';'#377eb8';'#ffffb3';'#a65628';'#ff7f00'};
% models = {'PA','Pt'}; G = 'PA';colorIDs = {'#a6d854';'#8da0cb'};
% models =  {'HT','PtHT','Sc','Yl','CHO'};G = 'HT';colorIDs = {'#4daf4a';'#377eb8';'#ffffb3';'#a65628';'#ff7f00'};

load(['3D_PA.mat']);
intMet = unique(mets(metTypes.AA));
colors = hex2rgb(colorIDs);
S = [];

for m = 1:length(models)
    load(['3D_',models{m},'.mat']);
    for e = 1:length(intMet)
        met = intMet{e};
        num = find(ismember(mets(group),met));
        if ~isempty(num)
            S(e,m) = groupSens(num);
            C(e,m) = groupCod(num);
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

figure('visible','on','rend','painters','pos',[10 10 700 500])
bar(S(:,1:end-1),'stacked')
xlim([0 length(intMet)+1])
set(gca, 'XTick', [1:length(intMet)])
set(gca, 'XTickLabel', textMet(indexS))
xlabel('Metabolites')
ylabel('Sensitivity')
legend(models,'Location','NorthWest')
colormap(colors)
saveas(gcf,['AAmets_S_ALL_',G],'svg')


figure('visible','on','rend','painters','pos',[10 10 700 500])
bar(C(:,1:end-1),'stacked')
xlim([0 length(intMet)+1])
set(gca, 'XTick', [1:length(intMet)])
set(gca, 'XTickLabel', textMet(indexC))
xlabel('Metabolites')
ylabel('Codependence')
legend(models,'Location','NorthWest')
colormap(colors)
saveas(gcf,['AAmets_C_ALL_',G],'svg')

%% ALL
clear
models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};

S = [];
for m = 1:length(models)
    load(['3D_',models{m},'.mat']);
    Mets{m} = metNums(group);
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
models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};
% models = {'Sc'}
S = [];
for m = 1:length(models)
    load(['3D_',models{m},'.mat']);
    Mets{m} = metNums(group);
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

%% Composition pie charts
models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};

for i = 1:length(models)
    load(['BOFsensitivity_',models{i}])
    [metTypes,metGroups,metNames] = findMetType(model,mets);
    F = compositionPieChart(model,mets,Stoich,metTypes,metGroups,ID);
    
end

%% Composition stacked bar
clear
models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};
colors = {'#a65628','#984ea3','#4daf4a','#ffff33','#ff7f00','#377eb8','#e41a1c'};
colors = hex2rgb(colors);

for m = 1:length(models)
    F = {};
    load(['BOFsensitivity_',models{m}])
    
    [metTypes,metGroups,metNames] = findMetType(model,mets);
    F = compositionPieChart(model,mets,Stoich,metTypes,metGroups,ID);
    
    compArr = table2cell(F);
    compMat{m} = cell2mat(compArr(:,2:end));
end

for n = 1:length(compMat{1}(:,1))
    for m = 1:length(models)
        typeMat(m,:) = compMat{m}(n,:);
    end
    STD(n,:) = std(typeMat,1);
    M(n,:) = mean(typeMat,1);
end
tM = [55;3.4+2.5+2.5;20.5+3.1;3.03;3.03;3.03;3.5+0.4]/100;
M = [tM,M];

figure
hold on
b = bar([0:length(compMat{1}(1,:))],M','stacked');
colormap(colors);

xpos = [1:length(compMat{1}(1,:))];
ypos = 0;
dx = zeros(1,length(compMat{1}(1,:)));

for n = 1:length(compMat{1}(:,1))
    dy = STD(n,:)/2;
    ypos = b(n).YData(2:end) + ypos;
    errorbar(xpos,ypos,dy,-dy,dx,-dx,...
        'LineStyle','none','LineWidth',1.1,'Color','k')
end

legend(compArr(:,1))
xlabel('Time')
ylabel('Mass fraction')
ylim([0 1.1])
box on

%% Composition group correlation
clear
models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};

for m = 1:length(models)
    F = {};
    load(['BOFsensitivity_',models{m}])
    
    [metTypes,metGroups,metNames] = findMetType(model,mets);
    F = compositionPieChart(model,mets,Stoich,metTypes,metGroups,ID);
    
    compArr = table2cell(F);
    compMat{m} = cell2mat(compArr(:,2:end));
end

corrMat = {};
for n = 1:length(compMat{1}(:,1))
    figure
    hold on
    x = [];
    y = [];
    for m = 1:length(models)
        t = length(compMat{1}(1,:));
        if all(compMat{m}(n,1:t))
            x = [x,1:t];
            y = [y,compMat{m}(n,1:t)/max(compMat{m}(n,1:t))];
        end
    end
    x = x';
    y = y';
    [fitobj,gof] = fit(x,y,'poly1');
    compArr{n,1}
    gof
    plot(fitobj,x,y)
%     STD(n,:) = std(typeMat,1);
%     M(n,:) = mean(typeMat,1);
end



%% C:N ratio 
clear
models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};
% models = {'PA'}
CNratio = {};

% fieldName = 'allMets';
fieldName = 'AA';
fitEq = 'exp1';

figure;
color = hex2rgb('#99d8c9');

for i = 1:length(models)
    load(['BOFsensitivity_',models{i}],'model')
    load(['3D_',models{i}])
    metNums = [];
    met = {};
    n = 0;
    subplot(2,4,i)
    hold on
    box on
    dgf = [];
    for m = 1:length(metTypes.(fieldName))
        temp = find(ismember(group,metTypes.(fieldName)(m)));
        if ~isempty(temp)
            n = n+1;
            metNums(n) = temp;
            met(n) = mets(metTypes.(fieldName)(m));
            metID = findMetIDs(model,met(n));
            [~,MW] = calculateFormula(model,metID,1);
            markerSize(n) = MW*100;
            formula = model.metFormulas{metID};
            Cnum = numAtomsOfElementInFormula(formula,'C');
            Nnum = numAtomsOfElementInFormula(formula,'N');
            
%             dgf(n) = DGf(m);
            
            CNratio{i}(n,1) = Nnum/Cnum;
%             CNratio{i}(n,1) = Nnum/Cnum/MW;           
        end
    end
    [xdata,index] = sort(CNratio{i});
    ydata = groupCost(metNums)';
    ydata = ydata(index);
    
    xdata = xdata/max(xdata);
    ydata = ydata/max(ydata);
    
    met = met(index);
    
    [coeffs,xfit,fdata,compData,finalError] = orthogonalFit(xdata,ydata,fitEq);

    d = 1.5*std(ydata);
    I = abs(finalError) >  d;
    
    % Fit without outliers
    [coeffs,xfit,fdata,compData,finalError] = orthogonalFit(xdata(I==0),ydata(I==0),fitEq);
    
%     dx = d*cos(pi/4);
    dx = 0;
    dy = d*sin(pi/4);
    
    plot(xdata,ydata,'o')
    plot(xfit,fdata)
    plot(xfit+dx,fdata+dy,'--b')
    plot(xfit-dx,fdata-dy,'--b')
    
    xlabel('$N:C$','Interpreter','latex')
    ylabel('Relative cost')
    xlim([0 1.05])
    ylim([0 1.05])
    
    metNames = cellfun(@(x) x(1:3),met,'UniformOutput',false);
    initial = cellfun(@(x) upper(x(1)),met,'UniformOutput',false);
    for j = 1:numel(initial)
        metNames{j}(1) = initial{j};
    end
    
    metNames(I==1)
    text(xdata,ydata,metNames)
    
end

%% Cost vs Gibbs free energy
clear
models = {'PA','HT','Pt','PtHT','CHO','Sc','Yl'};
% models = {'Yl','Sc'}
CNratio = {};

% fieldName = 'allMets';
fieldName = 'AA';
fitEq = 'poly1';

figure;
color = hex2rgb('#99d8c9');

aa_names = {'ala__L_c','arg__L_c','asn__L_c','asp__L_c','cys__L_c','gln__L_c','glu__L_c','gly_c','his__L_c','ile__L_c',...
            'leu__L_c','lys__L_c','met__L_c','phe__L_c','pro__L_c','ser__L_c','thr__L_c','trp__L_c','tyr__L_c','val__L_c'};
DGf = [-16.77;79.88;-43.53;-105.78;-9.31;-23.48;-82.27;-37.24;35;51.32;49.32;69.7;37.07;63.48;
        23.72;-50.2;-30.56;95.84;24.38;27.26];
    
% DGf = DGf - min(DGf);

for i = 1:length(models)
    metNames = [];
    metNums = [];
    load(['BOFsensitivity_',models{i}],'model')
    load(['3D_',models{i}])
    metNums = [];
    met = {};
    n = 0;
    subplot(2,4,i)
    hold on
    box on
    dgf = [];
    Nvec = [];
    for m = 1:length(metTypes.(fieldName))
        temp = find(ismember(group,metTypes.(fieldName)(m)));
        if ~isempty(temp)
            n = n+1;
            metNums(n) = temp;
            met(n) = mets(metTypes.(fieldName)(m));
            metID = findMetIDs(model,met(n));
            [~,MW] = calculateFormula(model,metID,1);
            markerSize(n) = MW*100;
            formula = model.metFormulas{metID};
            Cnum = numAtomsOfElementInFormula(formula,'C');
            Nnum = numAtomsOfElementInFormula(formula,'N');
            Nvec(metNums(n)) = Nnum;
            dgf(metNums(n)) = DGf(find(ismember(aa_names,met(n))));
            CNratio{i}(n,1) = Nnum/Cnum;
%             CNratio{i}(n,1) = Nnum/Cnum/MW;           
        end
    end
    [xdata,index] = sort(dgf(metNums)');
    ydata = groupCost(metNums)';
    ydata = ydata(index);
    
%     xdata = xdata/max(xdata);
%     ydata = ydata/max(ydata);
    
    met = met(index);
    
    [fitobj,gof] = fit(xdata,ydata,fitEq);
    rsq = gof.rsquare;
%     fitobj
    
%     [coeffs,xfit,fdata,compData,finalError] = orthogonalFit(xdata,ydata,fitEq);

%     d = 1.5*std(ydata);
%     I = abs(finalError) >  d;
    
    % Fit without outliers
%     [coeffs,xfit,fdata,compData,finalError] = orthogonalFit(xdata(I==0),ydata(I==0),fitEq);
%     dx = d*cos(pi/4);

%     dx = 0;
%     dy = d*sin(pi/4);

    coeffs = coeffvalues(fitobj);
%     models{i}
    ATPenergy = 1/coeffs(1);

    plot(fitobj,xdata,ydata,'o')
    legend('off')
%     plot(xfit+dx,fdata+dy,'--b')
%     plot(xfit-dx,fdata-dy,'--b')
    
    title(models{i},'FontName','Arial')
    xlabel('$\mathrm{\Delta G_{f}^{0}, kcal \cdot mole^{-1}}$','Interpreter','latex','FontName','Arial')
    ylabel('Biosynthetic cost','FontName','Arial')
    xlim([min(DGf)-5 max(DGf)]+5)
%     ylim([0 1.05])
    
    metNames = cellfun(@(x) x(1:3),met,'UniformOutput',false);
    initial = cellfun(@(x) upper(x(1)),met,'UniformOutput',false);
    for j = 1:numel(initial)
        metNames{j}(1) = initial{j};
    end
    
%     metNames(I==1)
    text(xdata+5,ydata,metNames,'FontName','Arial')
    text(0,max(ydata),['R2 = ',num2str(round(rsq,3))])
    
end

%% ATP energy vs. t
fieldName = 'AA';
group_costs = [];
models = {'PA','Pt','CHO','Sc','Yl'};
aa_names = {'ala__L_c','arg__L_c','asn__L_c','asp__L_c','cys__L_c','gln__L_c','glu__L_c','gly_c','his__L_c','ile__L_c',...
            'leu__L_c','lys__L_c','met__L_c','phe__L_c','pro__L_c','ser__L_c','thr__L_c','trp__L_c','tyr__L_c','val__L_c'};
DGf = [-16.77;79.88;-43.53;-105.78;-9.31;-23.48;-82.27;-37.24;35;51.32;49.32;69.7;37.07;63.48;
        23.72;-50.2;-30.56;95.84;24.38;27.26];
figure
for i = 1:length(models)
    i
    mets = [];
    Stoich = [];
    group_costs = [];
    group_mets = [];
    
    load(['BOFsensitivity_',models{i}])
    group_mets = [mets(metTypes.(fieldName))];
    bofMetIDs = findMetIDs(model,mets);
    group_costs = biosyntheticCost(model,group_mets,'atp_c -> adp_c + pi_c');   
    
    nums = find(ismember(aa_names,group_mets));
    xdata = DGf(nums');
    ydata = group_costs';
    fitobj = fit(xdata,ydata,'poly1');
    fitobj2 = fitlm(xdata,ydata)
    
    coeffs = coeffvalues(fitobj);
    ATPenergy = 1/coeffs(1)
    
    subplot(2,4,i)
    plot(fitobj,xdata,ydata,'o')
    text(xdata,ydata,group_mets)
%     title(models{i})
%     ylim([0 100])
end

%% Yield curves Biomass
modelsPA = {'PA','Pt','HT','PtHT'};
modelsHT = {'CHO','Yl','Sc'};
allModels = [modelsPA,modelsHT]

figure 
hold on
box on
Models = allModels;
for i = 1:(length(Models))
    Models{i}
    clear var yieldTc
    load(['BOFsensitivity_',Models{i}])
    if istable(yieldTc)
        yieldTc = table2cell(yieldTc);
    end
    pos = find(contains(yieldTc(:,1),'X'));
    
    allYc = cell2mat(yieldTc(:,2:end));
    corrCoeff = 1 + sum(allYc,1);
    Yt = cell2mat(yieldTc(pos,2:end));
    
    % Correct Yt for mass balance
    corrYt = Yt./corrCoeff

    %
    relYt = corrYt/max(corrYt);
    
    plot(relYt,'LineWidth',1)
    
    
end

legend(Models)

%% Yield curves X summary
modelsMicPA = {'PA','Pt'};
modelsMicHT = {'HT','PtHT'};
modelsHT = {'CHO','Yl','Sc'};
allModels = [modelsPA,modelsHT];

figure 
hold on
box on
Models = {modelsMicPA,modelsMicHT,modelsHT};
for i = 1:(length(Models))
    modelGroup = Models{i};
    for j = 1:length(modelGroup)
        clear var yieldTc
        groupYt = [];
        load(['BOFsensitivity_',modelGroup{j}])
        if istable(yieldTc)
            yieldTc = table2cell(yieldTc);
        end
        pos = find(contains(yieldTc(:,1),'X'));
        
        allYc = cell2mat(yieldTc(:,2:end));
        corrCoeff = 1 + sum(allYc,1);
        Yt = cell2mat(yieldTc(pos,2:end));
        
        % Correct Yt for mass balance
        corrYt = Yt./corrCoeff;
        
        groupYt = [groupYt;corrYt];
        %
        
    end
    plotYt = mean(groupYt,1);
    relYt = plotYt/max(plotYt);  
    plot(relYt,'LineWidth',1)
end

legend('MicPA','MicHT','HT')


%% Yield curves nitrogen
modelsPA = {'PA','Pt','HT','PtHT'};
modelsHT = {'CHO','Yl','Sc'};
allModels = [modelsPA,modelsHT]

figure 
hold on
box on
Models = allModels;
for i = 1:(length(Models))
    Models{i}
    clear var yieldTc
    load(['BOFsensitivity_',Models{i}])
    yieldTc = yieldT;
    if istable(yieldTc)
        yieldTc = table2cell(yieldTc);
    end
%     yieldTc
    pos = find(contains(yieldTc(:,1),'nh4_e'));
    if isempty(pos)
        pos = find(contains(yieldTc(:,1),'no3_e'));
    end
    Yt = abs(cell2mat(yieldTc(pos,2:end)));
    Yt = Yt/max(Yt);
    plot(Yt,'LineWidth',1)
    
    
end
legend(Models)

%% Yield curves N summary
modelsMicPA = {'PA','Pt'};
modelsMicHT = {'HT','PtHT'};
modelsHT = {'CHO','Yl','Sc'};
allModels = [modelsPA,modelsHT];

figure 
hold on
box on
Models = {modelsMicPA,modelsMicHT,modelsHT};
for i = 1:(length(Models))
    modelGroup = Models{i};
    for j = 1:length(modelGroup)
        clear var yieldTc
        groupYt = [];
        load(['BOFsensitivity_',modelGroup{j}])
        yieldTc = yieldT;
        if istable(yieldTc)
            yieldTc = table2cell(yieldTc);
        end
        pos = find(contains(yieldTc(:,1),'nh4_e'));
        if isempty(pos)
            pos = find(contains(yieldTc(:,1),'no3_e'));
        end
        
        Yt = abs(cell2mat(yieldTc(pos,2:end)));
        Yt = Yt/max(Yt);
        
        groupYt = [groupYt;Yt];
        %
        
    end
    plotYt = mean(groupYt,1);
    relYt = plotYt/max(plotYt);  
    plot(relYt,'LineWidth',1)
end

legend('MicPA','MicHT','HT')


%% Growth curves Biomass
modelsPA = {'PA','Pt','HT','PtHT'};
modelsHT = {'CHO','Yl','Sc'};
allModels = [modelsPA,modelsHT]

figure 
hold on
box on
Models = allModels;
uT = [];
for i = 1:(length(Models))

    load(['BOFsensitivity_',Models{i}])
    if istable(yieldTc)
        yieldTc = table2cell(yieldTc);
    end
    u = [];
    for i = 1:length(out)
        u(i) = out{i}.f;
    end
    u = u/max(u);
    plot(u,'LineWidth',1)
    
    uT = [uT;u];
    
end

legend(Models)


%% GP sampling
modelsPA = {'PA','Pt','HT','PtHT'};
modelsHT = {'CHO','Yl','Sc'};
Models = [modelsPA,modelsHT];

load(['BOFsensitivity_',Models{1}])
A = [length(Models),length(Stoich(1,:))];
SIZE = size(samplingResults);

for m = 1:A(1)
    model = [];
    load(['BOFsensitivity_',Models{m}])
    for i = 1:(length(Stoich(1,:)))
        if (m > SIZE(1)) || (i >= SIZE(2) && m == SIZE(1)) % To allow for resuming
            [m i]
            bof_id = find(model.c);
            model.S(findMetIDs(model,mets),bof_id) = -Stoich(:,i);
            [samplingResults{m,i}, mixedFrac(m,i)] = gpSampler(model);
            checkpoint = samplingResults{m,i};
        end
    end
end

save('samplingResults.mat','samplingResults', '-v7.3');

T = readtable('/home/jt/UCSD/BOFopt/pathways_manual.txt');
res = [];
for m = 1:length(Models)
    model = [];
    regulation = {};
    load(['BOFsensitivity_',Models{m}])
    for i = 1:length(model.rxns)
        x0 = samplingResults{m,1}.points(i,:);
        xf = samplingResults{m,end}.points(i,:);
        
        x0_mean = mean(x0);
        xf_mean = mean(xf);
        
        ttest = ttest(x0,xf);
        
        if ttest && xf_mean>x0_mean
            res(i,1) = 1;
        elseif  ttest && xf_mean<x0_mean
            res(i,1) = -1;
        else
            res(i,1) = 0;
        end
    end
    total = length(model.rxns);
    change = length(find(res));
    change/total
    regulation{m} = res;
end
x = samplingResults{1,1}.points(1300,:);
y = samplingResults{1,6}.points(1300,:);