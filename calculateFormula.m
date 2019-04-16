function [avMolarMass,molarMass,av,Num,formula] = calculateFormula(model,metIDs,weight)

%% Individual metabolite molecular weight

formulas = model.metFormulas(metIDs);
for m = 1:length(formulas)
    if isempty(model.metFormulas{metIDs(m)})
        warning(['No formula for ',model.mets{m}])
        molarMass(m,1) = 0;
    else
        Rnum(m,1) = numAtomsOfElementInFormula(formulas{m},'R');
        
        Cnum(m,1) = numAtomsOfElementInFormula(formulas{m},'C') + 16*Rnum(m,1);
        Hnum(m,1) = numAtomsOfElementInFormula(formulas{m},'H') + 32*Rnum(m,1);
        Onum(m,1) = numAtomsOfElementInFormula(formulas{m},'O') + 2*Rnum(m,1);
        Nnum(m,1)= numAtomsOfElementInFormula(formulas{m},'N');
        Pnum(m,1) = numAtomsOfElementInFormula(formulas{m},'P');
        Snum(m,1) = numAtomsOfElementInFormula(formulas{m},'S');
        
        molarMass(m,1) = (Cnum(m)*12+Hnum(m)+Onum(m)*16+Nnum(m)*14+Snum(m)*32+Pnum(m)*31)/1000;
    end
end

%% Average

avC = sum(Cnum.*weight)/sum(weight);
avH = sum(Hnum.*weight)/sum(weight);
avO = sum(Onum.*weight)/sum(weight);
avN = sum(Nnum.*weight)/sum(weight);
avP = sum(Pnum.*weight)/sum(weight);
avS = sum(Snum.*weight)/sum(weight);


av = [avC;avH;avO;avN;avP;avS];
av = av/avC;

Num = [Cnum,Hnum,Onum,Nnum,Pnum,Snum];
el = {'C','H','O','N','P','S'};
formula = '';
for f = 1:length(el)
    if av(f) == 1
        formula = [formula,el{f}];
    elseif av(f) < 5e-3
        formula = formula;
    else
        formula = [formula,el{f},num2str(round(av(f),2))];
    end
end

% disp(formula)
avMolarMass = sum(av.*[12;1;16;14;32;31])/1000; % g/mmol



end