
function [growth,parMetCoeff,dudp] = sensitivityAnalysis(model,met,mets,theta_0)

changeCobraSolver('gurobi','all',0);
fprintf('Calculating for #%4.0f (',m)
fprintf(char(parMet))
fprintf(')\n\n')

if range_Stoich(1) ~= range_Stoich(2)
    for b =1:itBOF
        for p = 1:itPoint
            model = changeRxnMets(model,mets,mets,bofRxn,theta_0(:,p,b));
            out=optimizeCbModel(model,'max','one');
            growth{m}(p,b)=out.f;
            if out.f == 0
                warning('Growth halted')
            end
        end
    end
    for b=1:itBOF
        for p = 1:itPoint-1
            dudp(p,b) = (growth(p+1,b)-growth(p,b))/(parMetCoeff(p+1,b)-parMetCoeff(p,b));
        end
    end
else
    growth = ones(itPoint,itBOF);
    parMetCoeff = ones(itPoint,itBOF);
    dudp = ones(itPoint-1,itBOF);
    fprintf(char(met))
    fprintf(' was excluded\n\n')
end

end