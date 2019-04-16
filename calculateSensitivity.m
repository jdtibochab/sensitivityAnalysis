%% Calculation of Biosynthetic Cost at constant BOF

function [S,Cod] = calculateSensitivity(model,mets,biomassReaction)
biomassReactionID = strmatch(biomassReaction,model.rxns);
bofMetIDs = find(model.S(:,biomassReactionID)~=0);
% mets = model.mets(bofMetIDs);
metIDs = find(ismember(model.mets,mets));
temp_model=model;
theta0 = temp_model.S(bofMetIDs,biomassReactionID);
change = [0.9:0.1:1.1]';
e=0;
r=0;
for m = 1:length(mets(:,1))
    fprintf('\nCalculating for ');fprintf(char(mets(m)));fprintf('\n\n')
    for b = 1:50
        metCoeffs(:,b) = theta0*change(1) + theta0.*rand(length(bofMetIDs),1).*(change(end)-change(1));
        for c = 1:length(change)
           metCoeff(c) = temp_model.S(metIDs(m),biomassReactionID)*change(c);
           metCoeffs(m,b) = metCoeff(c);
           temp_model.S(bofMetIDs,biomassReactionID)=metCoeffs(:,b);
           out = optimizeCbModel(temp_model,'max');
           if out.f>0
               r = r+1;
               metCoeff2(r) = metCoeff(c);
               response(r) = out.f;
           end
        end
            e = e+1;
            reg = polyfit(metCoeff2,response,1);
            Sb(e) = reg(1);
    end
    S(m) = mean(Sb); 
    Cod(m) = max(Sb)-min(Sb);
end