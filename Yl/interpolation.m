%% Interpolation
clear
[S,Mets] = xlsread('bofYl.xls');
S(1,:) = [];

[carbS,carbMets] = xlsread('carbYl.xls');
t = carbS(1,:)*24;

carbS(1,:) = [];
carbMets(1) = [];

t0 = [0 1 2 3 4 5 6 7]*24;
t = [12,39,87];

[intStoich,rsq] = interpolateStoich(carbS,t0,t);
