%% Parameter sensitivity and codependence analysis
clear
% Pt

fprintf('\n\n Calculating for Pt\n\n')

% load('ptData.mat')
% load('ptData2.mat')
load('ptData3.mat')
model_an = model_an2;

% [Stoech,Met]=xlsread('ptBOF.xls');
% Stoech=Stoech(2:end,2:end);
% Met = Met(2:354);
% [prodCapMet,notProdCapMet,Stoech2] = prodCap(modelPt,Met,Stoech);
[model_an,prodCapMet,notProdCapMet,Stoech2,theta_0] = metTest(modelPt,Met,Stoech);

Stoech = Stoech2;
Met = prodCapMet;

global Met_ID
Stoech = Stoech*-1;
%% Modification
biomassReactionIDs = [strmatch('BIOMASS_',model_an.rxns);strmatch('bof',model_an.rxns)];

model_an.ub(biomassReactionIDs(1:end-1)) = 0;

BOF = 'bof_c';

fprintf('\n \nRetrieving old and new stoichiometric coefficients\n\n')

model_analysis = model_an;

Met_ID=findMetIDs(model_analysis,Met);

fprintf('Generating optimization boundaries\n\n')
for i=1:length(Stoech(:,1))
    range_Stoech(i,:)=[min(Stoech(i,:)) max(Stoech(i,:))];
    init_Stoech(i,:)=[mean(Stoech(i,:))];
end

%%

it1 = 50; % BOF
it2 = 10; % Points

for k = 1:length(Met)
    parMet = Met(k);
    fprintf('Calculating for #%4.0f (',k)
    fprintf(char(parMet))
    fprintf(')\n\n')
    parMetID = findMetIDs(model_analysis,parMet);

    if range_Stoech(k,1) ~= range_Stoech(k,2)
        for j =1:it1
            theta_0(:,j)=range_Stoech(:,1)'+rand(1,length(range_Stoech(:,1))).*(range_Stoech(:,2)-range_Stoech(:,1))';
             for i = 1:it2
                parMetCoeff{k}(i,j) = range_Stoech(k,1)+(i/it2)*(range_Stoech(k,2)-range_Stoech(k,1));
                theta_0(k,j) = parMetCoeff{k}(i,j);
                model_analysis = changeRxnMets(model_analysis,model_analysis.mets(Met_ID),model_analysis.mets(Met_ID),BOF,theta_0(:,j));
                out=optimizeCbModel(changeObjective(model_analysis,BOF));
                growth{k}(i,j)=out.f;
                co2_uptake{k}(i,j)=out.full(strmatch('EX_glc__D_e',model_analysis.rxns));
             end
        end
        for j=1:length(theta_0(1,:))
           for i = 1:it2-1
               dudp{k}(i,j) = (growth{k}(i+1,j)-growth{k}(i,j))/(parMetCoeff{k}(i+1,j)-parMetCoeff{k}(i,j));
           end
        end
    else
        growth{k} = ones(it2,it1);
        parMetCoeff{k} = ones(it2,it1);
        dudp{k} = ones(it2-1,it1);
        co2_uptake{k} = ones(it2,it1);
        fprintf(char(parMet))
        fprintf(' was excluded\n\n')
    end

end
save('parSenCod_Pt.mat')

