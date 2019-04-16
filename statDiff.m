%% Statistical difference

function [p,t] = statDiff(M1,M2,S1,S2,N1,N2)

SE = sqrt(S1^2/N1 + S2^2/N2);

t = abs(M1-M2)/SE;

x1 = normrnd(M1,S1,N1,1);
x2 = normrnd(M2,S2,N2,1);

[h,p] = ttest2(x1,x2,'Vartype','unequal');

end