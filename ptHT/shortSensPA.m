clear
% load('ModelPA')
load('ptData3.mat')
model = model_an;
biomassReaction = find(model.c==1);


%%
Met = model.mets(find(model.S(:,strmatch(biomassReaction,model.rxns))));

AA_mets1 = strfind(Met,'__L');
AA_mets2 = strfind(Met,'gly_c');

j = 1;
for i = 1:length(AA_mets1)
   if AA_mets2{i} == 1
       AA_mets(j,1) = i;
       j = j+1;
   end
   if AA_mets1{i} == 4
       AA_mets(j,1) = i;
       j = j+1;
   end
end
%%

mets = AA_mets;
model = changeObjective(model,'bof_c');
%%
[Sens,Cod] = calculateSensitivity(model,Met(AA_mets),biomassReaction);

% Cost = biosyntheticCost(model,Met(group),0.19,'EX_glc-D_e'); % HT
% Cost = biosyntheticCost(model,Met(AA_mets),0.057,'EX_co2_e'); % Pt
Cost = biosyntheticCost(model,mets,0.025,'EX_co2(e)'); % PA

group = AA_mets;

%%
figure;
subplot(1,2,1)
plot(Cost,Sens,'*')
title('Sensitivity')
text(Cost+0.1,Sens,Met(group),'FontSize',8)
xlabel('A_{CO_2}')
ylabel('\Delta\mu/\Deltap')

subplot(1,2,2)
plot(Cost,Cod,'*')
title('Parameter codependence')
text(Cost+0.1,Cod,Met(group),'FontSize',8)
xlabel('A_{CO_2}')
ylabel('\Delta(\Delta\mu/\Deltap)')
x0=10;
y0=10;
width=1000;
height=800;
ylim([0 0.02])
set(gcf,'units','points','position',[x0,y0,width,height])