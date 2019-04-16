%% Parameter sensitivity and codependence analysis
clear
% AUTOTROPHY
load('ModelPAold')
load('ModelPA')

fprintf('\n\n Calculating for autotrophic growth\n\n')

model_an = [Model,ModelPA];

global Met_ID

for j=1:length(model_an)
    model=model_an{1,j};
    printRxnFormula(model,model.rxns(find(model.c==1)));
    solutionH{j}=optimizeCbModel(model);
    growth_initial(j) = solutionH{j}.f;
end
 % for each model you can extract the stoichiometry and associated
 % metabolites of the biomass equation
% Add the old data to the new one.
fprintf('\n \nRetrieving old and new stoichiometric coefficients\n\n')
for j=1:length(model_an)
    model=model_an{1,j};
    Srxn = full(model.S(:,find(model.c==1)));
        Srxn2{j} = Srxn(find(Srxn~=0));
    Met=model.mets(find(Srxn~=0));
        Met2{j} = Met;
end

for i=1:length(model_an)
    lengths(i)=length(Srxn2{i});
end
maxLength = max(lengths);

%% Vector length correction
for i=1:6
    for j=1:maxLength
       a=strcmp(Met2{i}(j),Met2{7}(j));
       if a == 0
           Met2{i}=[Met2{i}(1:j-1);Met2{7}(j);Met2{i}(j:end)];
           Srxn2{i}=[Srxn2{i}(1:j-1);0;Srxn2{i}(j:end)];
       end
    end   
end

Stoech = cell2mat(Srxn2);

model_analysis = model_an{1};

Met_ID=findMetIDs(model_analysis,Met);

fprintf('Generating optimization boundaries\n\n')
for i=1:length(Stoech(:,1))
    range_Stoech(i,:)=[min(Stoech(i,:)) max(Stoech(i,:))];
    init_Stoech(i,:)=[mean(Stoech(i,:))];
end

%%

it1 = 50; % BOF
it2 = 10; % Points

for k = 1:5 %length(Met)
    parMet = Met(k);
    fprintf('Calculating for #%4.0f (',k)
    fprintf(char(parMet))
    fprintf(')\n\n')
    parMetID = findMetIDs(model_analysis,parMet);

    if range_Stoech(k,1) ~= range_Stoech(k,2)
        for j =1:it1 % BOF
            theta_0(:,j)=range_Stoech(:,1)'+rand(1,length(range_Stoech(:,1))).*(range_Stoech(:,2)-range_Stoech(:,1))';
             for i = 1:it2 % Points
                parMetCoeff{k}(i,j) = range_Stoech(k,1)+(i/it2)*(range_Stoech(k,2)-range_Stoech(k,1));
                theta_0(k,j) = parMetCoeff{k}(i,j);
                model_analysis = changeRxnMets(model_analysis,model_analysis.mets(Met_ID),model_analysis.mets(Met_ID),'Biomass_Cvu_auto-',theta_0(:,j));
                out=optimizeCbModel(changeObjective(model_analysis,'Biomass_Cvu_auto-'));
                growth{k}(i,j)=out.f;
                co2_uptake{k}(i,j)=out.full(strmatch('EX_co2(e)',model_analysis.rxns));
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
save('parameterSenCod_PA.mat')


