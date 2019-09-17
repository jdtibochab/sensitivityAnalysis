
function [model,prodCapMets,notProdCapMets,Stoich,theta_0] = metTest(model,mets,Stoich)
bofID = find(model.c==1);

notBiomassMet = {'atp_c','adp_c','h2o_c','h_c','pi_c'};
notBiomassPos = [];
for i = 1:length(notBiomassMet)
    notBiomassPos = [notBiomassPos;strmatch(notBiomassMet{i},mets)];
end
biomassMet = mets;
biomassMet(notBiomassPos) = [];
biomassStoich = Stoich;
biomassStoich(notBiomassPos,:) = [];

Stoich = biomassStoich;
mets = biomassMet;
%% Check initial mass balance
for i = 1:length(Stoich(1,:))
    [~,compMolarMass] = calculateFormula(model,findMetIDs(model,mets),abs(Stoich(:,i)));
    testX = abs(sum(Stoich(:,i).*compMolarMass));
    if abs(testX-1) > 0.01
        warning('Initial S matrix does not meet mass balance. Mass balance will be forced upon the resulting matrix.')
    end
end

%% Identify metabolites to eliminate
metsInBOF = model.mets(find(model.S(:,bofID)));
eliminateMets = setdiff(metsInBOF,notBiomassMet);

elMetIDs = findMetIDs(model,eliminateMets);
model.S(elMetIDs,bofID) = 0;


%% Temporal model for calculations
temp_model = model;
temp_model.S(:,bofID) = 0;

%% Production capacity

MetIDs = findMetIDs(model,mets);
notProdCapMetNum = [];
  
model.S(MetIDs,bofID) = 0;

p=0;
n=0;
for j = 1:length(mets)
    metCoeff = -mean(Stoich(j,:));
    
    temp_model.S(MetIDs(j),bofID) = metCoeff;
    
    out = optimizeCbModel(temp_model,'max');
%     printRxnFormula(temp_model,temp_model.rxns(bofID));
%     disp(out.x(strmatch('DM_icit(h)',model.rxns)))
    
    if out.stat == 1 && out.f > 1e-4
       p = p+1;
       
       model.S(MetIDs(j),bofID) = metCoeff;
       
       prodCapMetIDs(p,1) = MetIDs(j);
       prodCapMetNum(p,1) = j;
       prodCapMets(p,1) = mets(j);
       disp([char(model.mets(MetIDs(j))),' included'])
    else
       n=n+1;
       
       model.S(MetIDs(j),bofID) = 0;
       temp_model.S(MetIDs(j),bofID) = 0;
       
       notProdCapMetIDs(n,1) = MetIDs(j);
       notProdCapMetNum(n,1) = j;
       notProdCapMets(n,1) = mets(j);

       disp([char(model.mets(MetIDs(j))),' discarded'])
    end
end

for i = 1 : length(notProdCapMetNum)
    j = length(notProdCapMetNum)-i+1;
    Stoich(notProdCapMetNum(j),:) = [];
end

prodCapMets = mets(prodCapMetNum);
notProdCapMets = mets(notProdCapMetNum);

%% Negative corrections
Stoich(Stoich<0) = 0;

%% Biomass check for 1g
for i = 1:length(Stoich(1,:))
    [molarX,compMolarMass,av,Num,formula] = calculateFormula(model,findMetIDs(model,prodCapMets),abs(Stoich(:,i)));
    testX = abs(sum(Stoich(:,i).*compMolarMass));
    Stoich(:,i) = Stoich(:,i)/testX;
end

%%
theta_0 = model.S(find(model.S(:,bofID)<0),bofID);
printRxnFormula(model,model.rxns(bofID));
end
