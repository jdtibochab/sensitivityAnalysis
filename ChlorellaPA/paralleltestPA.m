

%% Parameter sensitivity and codependence analysis
clear
tic

load('ModelPAold')
load('ModelPA')

fprintf('\n\n Calculating for autotrophic growth\n\n')

model_an = [Model,ModelPA];
for m = 1:length(model_an)
%    model_an{m} = starchCorr(model_an{m},'starch300_h','DM_o2D(u)',0.0355);
   out = optimizeCbModel(model_an{m});
   out.f
end
%%
solutionH = cell(1,length(model_an));
growth_initial = zeros(length(model_an));
for b=1:length(model_an)
    model=model_an{1,b};
    printRxnFormula(model,model.rxns(find(model.c==1)));
    solutionH{b}=optimizeCbModel(model);
    growth_initial(b) = solutionH{b}.f;
end
    
%%
fprintf('\n \nRetrieving old and new stoichiometric coefficients\n\n')

Srxn2 = cell(1,length(model_an));
Met2 = cell(1,length(model_an));
for b=1:length(model_an)
    model=model_an{1,b};
    Srxn = full(model.S(:,find(model.c==1)));
        Srxn2{b} = Srxn(find(Srxn~=0));
    Met=model.mets(find(Srxn~=0));
        Met2{b} = Met;
end

lengths = zeros(length(model_an));
for p=1:length(model_an)
    lengths(p)=length(Srxn2{p});
end
maxLength = max(lengths);

%% Vector length correction
for p=1:6
    for b=1:maxLength
       a=strcmp(Met2{p}(b),Met2{7}(b));
       if a == 0
           Met2{p}=[Met2{p}(1:b-1);Met2{7}(b);Met2{p}(b:end)];
           Srxn2{p}=[Srxn2{p}(1:b-1);0;Srxn2{p}(b:end)];
       end
    end   
end

Stoich = cell2mat(Srxn2);

%%
t0 = [3 4 5 6 7 8];
t = [4 4.5 5 5.5 6 6.5];
oldStoich = Stoich(:,1:6);
newStoich = Stoich(:,7:end);
[oldStoich,rsq] = interpolateStoich(oldStoich,t0,t);
[newStoich,rsq] = interpolateStoich(newStoich,t0,t);
Stoich = [oldStoich,newStoich]; 
% Stoich = (oldStoich+newStoich)/2;

%% Biomass composition change
N = 1;
mets = Met2{1};
for i = 1:length(Stoich(1,:))
    molarCoeffs = abs(Stoich(:,i));
    molarCoeffs(molarCoeffs>1) = 0;
    metIDs = findMetIDs(model_an{1},mets);
    [~,molarMass] = calculateFormula(model,metIDs,1);
    massCoeffs = molarCoeffs.*molarMass;
    totalMass(i) = sum(massCoeffs);
    massFrac{i} = massCoeffs/totalMass(i);
end
[metTypes,mets] = findMetType(model,mets);

metGroups.CB = {'CB'};
metGroups.rest = {'rest'};
metGroups.Nuc = {'Nuc'};
metGroups.AA = {'AA'};
metGroups.TAG = {'TAG'};
metGroups.PL = {'Pe','Pg','Pail','Pchol'};
metGroups.OL = {'As','Dgts','Dgdg','Sqdg','Mgdg'};
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
%     subplot(N,length(Stoich(1,:))/N,i)
%     pie(fraction(:,i),metFieldNames)
end
%%
% otherMets = {'atp_c';'adp_c';'h_c';'pi_c'};
% changes = 0.99 + (1.01-0.99)*rand(1,length(Stoech(1,:)));
% otherMetNum = find(ismember(Met,otherMets));
% 
% for i = 1:length(otherMetNum)
%    Stoech(otherMetNum(i),:) =  Stoech(otherMetNum(i),:).*changes;
% end

model_analysis = model_an{1};

Met_ID=findMetIDs(model_analysis,Met);

fprintf('Generating optimization boundaries\n\n')

range_Stoich = zeros(length(Met),2);
for p=1:length(Met)
    range_Stoech(p,:)=[min(Stoich(p,:)) max(Stoich(p,:))];
end



%%
clear out
% rN = [6.5,(6.5+2.1)/2,((6.5+2.1)/2+2.1)/2, 2.1, 0.5, 0.25]/6.5*-0.2114;
model_temp = model_an{1};


biomassMets = Met;
biomassStoich = Stoich;
eliminateMets = {'h_c','pi_c','atp_c','h2o_c','adp_c'};

for i = 1:length(eliminateMets)
    elPos = find(ismember(biomassMets,eliminateMets{i}));
    biomassMets(elPos) = [];
    biomassStoich(elPos,:) = [];
end

for i = 1:length(biomassStoich(1,:))
    bofMetIDs = findMetIDs(model_temp,biomassMets);
%     model_temp.lb(strmatch('EX_no3(e)',model_temp.rxns)) = rN(i);
    model_temp.S(bofMetIDs,find(model_temp.c==1)) = biomassStoich(:,i);
    modelsPA{i} = model_temp;
    out{i} = optimizeCbModel(model_temp,'max','one');
    u(i) = out{i}.f;
    
    [molarX,compMolarMass,av,Num,formula] = calculateFormula(model_temp,findMetIDs(model_temp,biomassMets),abs(biomassStoich(:,i)));
    testX = abs(sum(biomassStoich(:,i).*compMolarMass));
    
    rS = abs(out{i}.x(strmatch('EX_co2(e)',model_temp.rxns)))*1;
    rX = abs(out{i}.f)/molarX*av(1);
    rC = 0;
    rN = abs(out{i}.x(strmatch('EX_no3(e)',model_temp.rxns)));
    ro = 0;
    
    Ysa(i) = ro/(rS);
    Ysx(i) = rX/(rS);
    Ysc(i) = rC/(rS);
    test(i) = Ysa(i)+Ysx(i)+Ysc(i);
    
    Ynx(i) = rX/rN; % C-molX/molN;
    
end
R = [Ysa;Ysx;Ysc;Ynx];

%%
CSVfile = 'fluxDist_PA';
writeFluxCsv(model_temp,out,CSVfile);

%% Variations
clear deltaU deltaZ rxnID

rxnID = [1:length(model.rxns)]'*ones(1,length(u));
el = [];
for i = 2:length(Stoich(1,:))
    deltaU(:,i) = abs((out{i}.x - out{1}.x));
end
deltaU(:,1) = [];
deltaU(abs(deltaU)>100) = 0;

for i = 1:length(deltaU(:,1))
    if sum(abs(deltaU(i,:))) == 0
        el(end+1) = i;
    end
end
rxnID(el,:) = [];
deltaU(el,:) = [];
rxnID = rxnID(:,1);

deltaZ = zscore(deltaU,0,2);

%%
it1 = 100; % BOF
it2 = 10; % Points

theta_0 = cell(length(Met),1);
parMetCoeff = cell(length(Met),1);
for m=1:length(Met)
    if range_Stoech(m,1) ~= range_Stoech(m,2)
        for b=1:it1 % BOF
             theta_0{m}(:,:,b)=(range_Stoech(:,1)+rand(length(range_Stoech(:,1)),1).*(range_Stoech(:,2)-range_Stoech(:,1)))*ones(1,it2);
             for p=1:it2 % Point
                 parMetCoeff{m}(p,b) = range_Stoech(m,1)+(p/it2)*(range_Stoech(m,2)-range_Stoech(m,1));
                 theta_0{m}(m,p,b) = parMetCoeff{m}(p,b);
             end
        end
    end
end
growth = cell(length(Met),1);

% parpool(12)
parfor m = 1:length(Met)
    changeCobraSolver('gurobi','all',0);
    model_analysis = model_an{1};
    parMet = Met(m);
    fprintf('Calculating for #%4.0f (',m)
    fprintf(char(parMet))
    fprintf(')\n\n')    
    parMetID = findMetIDs(model_analysis,parMet);
    if range_Stoech(m,1) ~= range_Stoech(m,2)
        for b =1:it1
            for p = 1:it2
                model_analysis = changeRxnMets(model_analysis,model_analysis.mets(Met_ID),model_analysis.mets(Met_ID),'Biomass_Cvu_auto-',theta_0{m}(:,p,b));
                out=optimizeCbModel(changeObjective(model_analysis,'Biomass_Cvu_auto-'),'max','one');
                growth{m}(p,b)=out.f;
                co2_uptake{m}(p,b)=out.x(strmatch('EX_co2(e)',model_analysis.rxns));
            end
        end
        for b=1:it1
           for p = 1:it2-1
               dudp{m}(p,b) = (growth{m}(p+1,b)-growth{m}(p,b))/(parMetCoeff{m}(p+1,b)-parMetCoeff{m}(p,b));
           end
        end
    else
        growth{m} = ones(it2,it1);
        parMetCoeff{m} = ones(it2,it1);
        dudp{m} = ones(it2-1,it1);
        co2_uptake{m} = ones(it2,it1);
        fprintf(char(parMet))
        fprintf(' was excluded\n\n')
    end

end
time=toc/60;
fprintf('\n\n Calculation time = %4.2f min \n\n',time)
save('BOFsensitivity_PA.mat')