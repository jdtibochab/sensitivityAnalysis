%
%Model fixed ready for TestSolutionChlorella
%No loopps
%
%
%
clear all
initCobraToolbox    

load('C:\Users\Cristal\Documents\MATLAB\Chlorella1\manual_curation\Workspaces\iCZ842_ShadowPrices.mat')
DW = 13*10^(-12); % avg. dry weight of log phase chlamy cell = 13 pg (This study)
CPerStarch300 = 1800; % derived from starch300 chemical formula
starchPAdD=0.042; % Measured in dark g/gDW d
starchPAd=0.113; %Measured in light g/gDW d
starchDegAerLight = starchPAd/(CPerStarch300*24);% estimated fron the rates calculated with experimental data mmol/g/h
starchDegAnLight = (3/2)*starchDegAerLight; % approx. SS rate of aerobic starch degradation in light 3/2 of anaerobic rate in green algae (Nakamura and Miyachi, 1982)
starchDegAnDark = starchPAdD/(CPerStarch300*24); % estimated fron the rates calculated with experimental data mmol/g/h
starchDegAerDark = (3/2)*starchDegAnDark; % approx. SS rate of aerobic starch degradation in light 3/2 of anaerobic rate in green algae (Nakamura and Miyachi, 1982)
dimensionalConversion = 3.836473679; % from emitted microE/m^2/s to incident mmol/gDW/hr
effectiveConversion = 0.06; % from incident mmol/gDw/hr to effective mmol/gDw/hr [Hiramaya 1995]
% model=cvuI24;

model = changeRxnBounds(model,{'MALOAAtm'},0,'b');%PRevent loops
model = changeRxnBounds(model,{'OAACITtm'},0,'l');%PRevent loops
model = changeRxnBounds(model,{'OAAAKGtm'},0,'l');%PRevent loops
model = changeRxnBounds(model,{'FADth'},0,'u');%PRevent loops
model = changeRxnBounds(model,{'PYRtm'},0,'u');%PRevent loops

save('iCZ842_28Sep','model');

%Photoautotrophy
modelLna = model;
% The single PRISM reaction being used has to be commented-out below.
modelLna = changeRxnBounds(modelLna,{...
%    'PRISM_solar_litho',...
    'PRISM_design_growth',...
    'PRISM_solar_exo',...
    'PRISM_fluorescent_warm_18W',...
    'PRISM_incandescent_60W',...
    'PRISM_fluorescent_cool_215W',...
    'PRISM_metal_halide',...
    'PRISM_high_pressure_sodium',...
    'PRISM_growth_room',...
    'PRISM_white_LED',...
    'PRISM_red_LED_array_653nm',...
    'PRISM_red_LED_674nm',...
},0,'b');
modelLna = changeRxnBounds(modelLna,{'EX_o2(e)'},-1000,'l');
modelLna = changeRxnBounds(modelLna,{'EX_co2(e)'},-13.54,'l');
modelLna = changeRxnBounds(modelLna,{'EX_ac(e)'},0,'l');
modelLna = changeRxnBounds(modelLna,{'EX_no3(e)'},-10,'l');
modelLna = changeRxnBounds(modelLna,{'EX_nh4(e)'},0,'l');
modelLna = changeRxnBounds(modelLna,{'EX_starch(h)'},0,'b');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGRA',starchDegAerLight,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGR2A',0,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGRB',starchDegAerLight,'u');
modelLna = changeRxnBounds(modelLna,'STARCH300DEGR2B',0,'u');
modelLna = changeRxnBounds(modelLna,{'PCHLDR'},0,'b'); % Crystal Structure of the Nitrogenase-like Dark Operative Protochlorophyllide Oxidoreductase Catalytic Complex [Brocker 2010]
modelLna = changeRxnBounds(modelLna,{'G6PADHh','G6PBDHh'},0,'b'); % Purication of chloroplast G6PDH had not been successful, partly because this isoform is inactivated by light because of redox modulation by thioredoxin, and partly because it easily aggregates unspecically during purication in C. vulgaris [Honjoh, 2003].
modelLna = changeRxnBounds(modelLna,{'FBAh'},0,'b'); % light inactivates FBAh [Willard and Gibbs, 1968]
modelLna = changeRxnBounds(modelLna,{'H2Oth'},0,'u'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark (Packer 1970)
modelLna = changeRxnBounds(modelLna,{'Biomass_Cvu_mixo-','Biomass_Cvu_hetero-'},0,'b');
modelLna = changeObjective(modelLna,'Biomass_Cvu_auto-');
%Solution
solutionLna = optimizeCbModel(modelLna);
solutionLna.f

%%
%%% light, aerobic, CO2, w/glucose, biomass objective
modelLMix = model;
% The single PRISM reaction being used has to be commented-out below.
modelLMix = changeRxnBounds(modelLMix,{...
%    'PRISM_solar_litho',...
    'PRISM_solar_exo',...
    'PRISM_incandescent_60W',...
    'PRISM_fluorescent_cool_215W',...
    'PRISM_metal_halide',...
    'PRISM_high_pressure_sodium',...
    'PRISM_growth_room',...
    'PRISM_white_LED',...
    'PRISM_red_LED_array_653nm',...
    'PRISM_red_LED_674nm',...
    'PRISM_fluorescent_warm_18W',...
    'PRISM_design_growth',...    
},0,'b');
modelLMix = changeRxnBounds(modelLMix,{'EX_o2(e)',},-1000,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_glc-A(e)'},-0.3025,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_ac(e)'},-0,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_starch(h)'},0,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_no3(e)'},-100,'l');
modelLMix = changeRxnBounds(modelLMix,{'EX_nh4(e)'},0,'l');
modelLMix = changeRxnBounds(modelLMix,'EX_co2(e)',-13.6,'l');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGRA',starchDegAerLight,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGR2A',0,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGRB',starchDegAerLight,'u');
modelLMix = changeRxnBounds(modelLMix,'STARCH300DEGR2B',0,'u');
modelLMix = changeRxnBounds(modelLMix,{'PCHLDR'},0,'b'); % Crystal Structure of the Nitrogenase-like Dark Operative Protochlorophyllide Oxidoreductase Catalytic Complex [Brocker 2010]
modelLMix = changeRxnBounds(modelLMix,{'G6PADHh','G6PBDHh'},0,'b'); % Purication of chloroplast G6PDH had not been successful, partly because this isoform is inactivated by light because of redox modulation by thioredoxin, and partly because it easily aggregates unspecically during purication in C. vulgaris [Honjoh, 2003].
modelLMix = changeRxnBounds(modelLMix,{'FBAh'},0,'b'); % light inactivates FBAh [Willard and Gibbs, 1968]
modelLMix = changeRxnBounds(modelLMix,{'H2Oth'},0,'u'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark [Packer 1970]
modelLMix = changeRxnBounds(modelLMix,{'Biomass_Cvu_auto-','Biomass_Cvu_hetero-'},0,'b');
modelLMix = changeObjective(modelLMix,'Biomass_Cvu_mixo-');
solutionLMix = optimizeCbModel(modelLMix);
solutionLMix.f

%%
%%% dark, aerobic, w/ glucose, biomass objective
modelDa = model;
modelDa = changeRxnBounds(modelDa,'EX_photonVis(e)',0,'l');
modelDa = changeRxnBounds(modelDa,{'PRISM_solar_litho','PRISM_solar_exo','PRISM_incandescent_60W','PRISM_fluorescent_cool_215W','PRISM_metal_halide','PRISM_high_pressure_sodium','PRISM_growth_room','PRISM_white_LED','PRISM_red_LED_array_653nm','PRISM_red_LED_674nm','PRISM_fluorescent_warm_18W','PRISM_design_growth'},0,'b');
modelDa = changeRxnBounds(modelDa,{'EX_o2(e)'},-10,'l');
modelDa = changeRxnBounds(modelDa,'EX_hco3(e)',-0,'l');
modelDa = changeRxnBounds(modelDa,'EX_ac(e)',-0,'l');
modelDa = changeRxnBounds(modelDa,'EX_no3(e)',-10,'l');
modelDa = changeRxnBounds(modelDa,'EX_nh4(e)',0,'l');
modelDa = changeRxnBounds(modelDa,'EX_h(e)',-10,'l');
modelDa = changeRxnBounds(modelDa,'EX_for(e)',0,'b');
modelDa = changeRxnBounds(modelDa,'EX_succ(e)',0,'b');
modelDa = changeRxnBounds(modelDa,'EX_glc-A(e)',-0.3025,'l');
modelDa = changeRxnBounds(modelDa,'STARCH300DEGRA',0,'u');
modelDa = changeRxnBounds(modelDa,'STARCH300DEGR2A',starchDegAerDark,'u');
modelDa = changeRxnBounds(modelDa,'STARCH300DEGRB',0,'u');
modelDa = changeRxnBounds(modelDa,'STARCH300DEGR2B',starchDegAerDark,'u');
modelDa = changeRxnBounds(modelDa,{'ATPSh'},0,'b'); % Inactive in dark It is suggested that there exists in dark-adapted algae a permanent proton gradient which stimulates the charge recombination process. This proton gradient results from the hydrolysis of a pool of ATP by membranebound ATPases and collapses after the addition of TNBT. The long lifetime of this proton gradient (several hours) indicates that the ATP probably comes from the mitochondria [Joliot and Joliot, 1980].
modelDa = changeRxnBounds(modelDa,{'GAPDH(nadp)hi'},0,'b'); % Only active in light [Buchanan 1980]
modelDa = changeRxnBounds(modelDa,{'MDH(nadp)hi','MDHC(nadp)hr'},0,'b'); % inactive in dark [Buchanan 1980]
modelDa = changeRxnBounds(modelDa,{'PPDKh'},0,'b'); % inactive in dark [Buchanan 1980]
modelDa = changeRxnBounds(modelDa,{'PRUK'},0,'b'); % inactive in dark [Buchanan 1980; Villarejo 1995]
modelDa = changeRxnBounds(modelDa,{'RBPCh','RBCh'},0,'b'); % That Rubisco was first identified as a light-activated enzyme in early labeling studies with air-grown Chlorella cells provides evidence in support of this conclusion [Pedersen et al., 1966].
modelDa = changeRxnBounds(modelDa,{'SBP'},0,'b'); % inactive in dark [Buchanan 1980, Koziol 2013] 
modelDa = changeRxnBounds(modelDa,{'H2Oth'},0,'l'); % there is a high h2o requirement in [h]; however, experiments show that h2o in general goes from [h] to [c] in light and from [c] to [h] in dark [Packer 1970]
modelDa = changeRxnBounds(modelDa,{'Biomass_Cvu_auto-','Biomass_Cvu_mixo-'},0,'b');
modelDa = changeObjective(modelDa,'Biomass_Cvu_hetero-');
solutionDa = optimizeCbModel(modelDa);
solutionDa.f



spot={'MALOAAtm';'ASPAT';'OAAAKGtm';'OAACITtm';'CSm';'PYRtm'};
results={};
for i=1:numel(spot)
    a=strmatch(spot{i},model.rxns,'exact');
    if i==1
    results{1,1}=solutionLna.x(a);
    results{2,1}=solutionLMix.x(a);
    results{3,1}=solutionDa.x(a);
    elseif i==2
    results{1,2}=solutionLna.x(a);
    results{2,2}=solutionLMix.x(a);
    results{3,2}=solutionDa.x(a);
    elseif i==3
    results{1,3}=solutionLna.x(a);
    results{2,3}=solutionLMix.x(a);
    results{3,3}=solutionDa.x(a);
    elseif i==4
    results{1,4}=solutionLna.x(a);
    results{2,4}=solutionLMix.x(a);
    results{3,4}=solutionDa.x(a);
    elseif i==5
    results{1,5}=solutionLna.x(a);
    results{2,5}=solutionLMix.x(a);
    results{3,5}=solutionDa.x(a);
    elseif i==6
    results{1,6}=solutionLna.x(a);
    results{2,6}=solutionLMix.x(a);
    results{3,6}=solutionDa.x(a);
    end 
end


%%
%Editing iCZ842
%
%
%




