%% BOF sensitivity analysis
% Computes the sensitivity of a certain model to the change of metabolite
% concentration in the biomass. This code calculates the sensitivity as the
% change of growth rate as a result of the change in the stoichiometric
% coefficient of each metabolite in [mets]. The calculation includes the
% estimation of the codependence as the change of sensitivity depending on
% the composition of other metabolites.

% This calculation might take long, typically from several minutes to hours.

function [growth,parMetCoeff,dudp] = BOFsensitivity(model,mets,calcMets,Stoich,itBOF,itPoint,N)

Stoich = -Stoich;
metIDs = findMetIDs(model,mets);
bofID = find(model.c==1);
bofRxn = model.rxns(bofID);

for m = 1:length(calcMets)
    calcMetNum(m) = find(ismember(mets,calcMets{m}));
end

range_Stoich = zeros(length(mets),2);
for p = 1:length(mets)
    range_Stoich(p,:)=[min(Stoich(p,:)) max(Stoich(p,:))];
end

%% Sensitivity analysis

theta_0 = cell(length(mets),1);
parMetCoeff = cell(length(mets),1);
for m = 1:length(mets)
    if range_Stoich(m,1) ~= range_Stoich(m,2)
        for b=1:itBOF % BOF
             theta_0{m}(:,:,b)=(range_Stoich(:,1)+rand(length(range_Stoich(:,1)),1).*(range_Stoich(:,2)-range_Stoich(:,1)))*ones(1,itPoint);
             for p=1:itPoint % Point
                 parMetCoeff{m}(p,b) = range_Stoich(m,1)+(p/itPoint)*(range_Stoich(m,2)-range_Stoich(m,1));
                 theta_0{m}(m,p,b) = parMetCoeff{m}(p,b);
             end
        end
    end
end
growth = cell(length(mets),1);
parModel = model;


delete(gcp('nocreate'))
parpool(N)

parfor m = 1:length(mets)
    model = parModel;
    changeCobraSolver('gurobi','all',0);
    parMet = mets(m);
    parMetID = findMetIDs(model,parMet);
    if range_Stoich(m,1) ~= range_Stoich(m,2) && ismember(m,calcMetNum)
        disp(['Calculating for #',num2str(m),' (',char(parMet),')'])
        for b =1:itBOF
            for p = 1:itPoint
                model = changeRxnMets(model,mets,mets,bofRxn,theta_0{m}(:,p,b));
                out=optimizeCbModel(model,'max','one');
%                 out=optimizeCbModel(model);
                growth{m}(p,b)=out.f;
                if out.f == 0
                    warning('Growth halted')
                end
            end
        end
        for b = 1:itBOF
            x = parMetCoeff{m}(:,b);
            y = growth{m}(:,b);
            
            poly = polyfit(x,y,1);
            
            dudp{m}(b) = poly(1);
        end
    else
        dudp{m} = zeros(1,itBOF);
        growth{m} = ones(itPoint,itBOF);
        parMetCoeff{m} = ones(itPoint,itBOF);
    end
delete(gcp('nocreate'))
end