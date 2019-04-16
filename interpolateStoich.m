function [intStoich,rsq] = interpolateStoich(stoich,t0,t)
maxDeg = length(t0)-1;
maxDeg = min([maxDeg,3]);
x = t0;
for m = 1:length(stoich(:,1))
    y = stoich(m,:);
    for i = 1:maxDeg
        poly = ['poly',num2str(i)];
        [polyfit{i},gof] = fit(x',y',poly);
        R2(i) = gof.adjrsquare;
    end
    R2(isnan(R2)) = 0;
    [rsq(m),deg] = max(R2);
    intStoich(m,:) = feval(polyfit{deg},t);
end
end