clear
% Pt
tic
fprintf('\n\n Calculating for Pt TEST1\n\n')

% load('ptData.mat')
% load('ptData3.mat')

% load('modelPt.mat')

load('ptData_norm.mat')
biomassReactionIDs = [strmatch('BIOMASS_',model_an.rxns);strmatch('bof',model_an.rxns)];
model_an.lb(biomassReactionIDs(1:end-1)) = 0;
model_an.ub(biomassReactionIDs(1:end-1)) = 0;

% 
[Stoich,Met]=xlsread('ptBOFnorm.xls');
modelPt = changeObjective(modelPt,'bof_c');

%% HETEROTROPHY
% Activate 
% modelPt = addExchangeRxn(modelPt,'glyc_e');
modelPt = changeRxnBounds(modelPt,'EX_glc__D_e',-1000,'l');
% modelPt = addReaction(modelPt,'GLYt','reactionFormula','glyc_e -> glyc_c');
% modelPt.csense = [modelPt.csense,'E'];

% Carbon exchange
modelPt = changeRxnBounds(modelPt,'EX_co2_e',0,'l');
modelPt = changeRxnBounds(modelPt,'EX_co2_e',1000,'u');

% modelPt = changeRxnBounds(modelPt,'EX_lac__D_e',-1000,'l');

optimizeCbModel(modelPt);
[a,b] = exchangeSingleModel(modelPt);

%%

% [model_an,prodCapMet,notProdCapMet,Stoech2] = prodCap(modelPt,Met,Stoech);

[model_an,prodCapMet,notProdCapMet,Stoech2,theta_0] = metTest(modelPt,Met,Stoich);
% % 
Stoich = Stoech2;
Met = prodCapMet;

global Met_ID
Stoich = Stoich*-1;

t0 = [1	1.5	2	3	5	12	17];
t = [4 4.5 5 5.5 6 6.5];
[intStoich,rsq] = interpolateStoich(Stoich,t0,t);
Stoich = intStoich;

%%
% exchange = 'EX_no3_e';
% exchange = 'EX_o2_e';
exchange = 'EX_photon_e';
exchangeC = 'EX_glc__D_e';

model_an.lb(strmatch(exchangeC,model_an.rxns)) = -13.54;


[model_an,reqExchange] = fixGrowth(model_an,exchange,0.057);

model_an.lb(strmatch(exchange,model_an.rxns)) = reqExchange;

model_an.lb(find(model_an.c==1)) = modelPt.lb(find(model_an.c==1));
model_an.ub(find(model_an.c==1)) = modelPt.ub(find(model_an.c==1));

[a,b] = exchangeSingleModel(model_an);
out = optimizeCbModel(model_an)


rN = [6.5,(6.5+2.1)/2,((6.5+2.1)/2+2.1)/2, 2.1, 0.5, 0.25, 0]/6.5*-1.76;
model_temp = model_an;
for i = 1:length(Stoich(1,:))
    bofMetIDs = findMetIDs(model_temp,Met);
    model_temp.lb(strmatch('EX_no3_e',model_temp.rxns)) = rN(i);
    model_temp.S(bofMetIDs,find(model_temp.c==1)) = Stoich(:,i);
    modelsPt{i} = model_temp;
    out = optimizeCbModel(model_temp);
    u(i) = out.f;
end

%% Modification
biomassReactionIDs = [strmatch('BIOMASS_',model_an.rxns);strmatch('bof',model_an.rxns)];

model_an.ub(biomassReactionIDs(1:end-1)) = 0;

BOF = model_an.rxns(find(model_an.c==1));
% model = changeObjective(model,BOF);

fprintf('\n \nRetrieving old and new stoichiometric coefficients\n\n')

model_analysis = model_an;

Met_ID=findMetIDs(model_analysis,Met);

fprintf('Generating optimization boundaries\n\n')
for i=1:length(Stoich(:,1))
    range_Stoech(i,:)=[min(Stoich(i,:)) max(Stoich(i,:))];
    init_Stoech(i,:)=[mean(Stoich(i,:))];
end

%%

it1 = 50; % BOF
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

parpool(12)
parfor m = 1:length(Met)
    changeCobraSolver('gurobi','all',0);
    model_analysis = model_an;
    parMet = Met(m);
    fprintf('Calculating for #%4.0f (',m)
    fprintf(char(parMet))
    fprintf(')\n\n')    
    parMetID = findMetIDs(model_analysis,parMet);
    if range_Stoech(m,1) ~= range_Stoech(m,2)
        for b =1:it1
            for p = 1:it2
                model_analysis = changeRxnMets(model_analysis,model_analysis.mets(Met_ID),model_analysis.mets(Met_ID),BOF,theta_0{m}(:,p,b));
                out=optimizeCbModel(changeObjective(model_analysis,BOF));
                growth{m}(p,b)=out.f;
%                 co2_uptake{m}(p,b)=out.full(strmatch('EX_glyc_e',model_analysis.rxns));
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
%         co2_uptake{m} = ones(it2,it1);
        fprintf(char(parMet))
        fprintf(' was excluded\n\n')
    end

end
time=toc/60;
fprintf('\n\n Calculation time = %4.2f min \n\n',time)
save('paralleltestPtHT_interpol.mat')