function J= fcoutBiomass(theta_min)
global Met_ID model_opt G_data

    [model_opt] = changeRxnMets(model_opt,model_opt.mets(Met_ID),model_opt.mets(Met_ID),'Biomass_Cvu_hetero-',theta_min);

    out=optimizeCbModel(changeObjective(model_opt,'Biomass_Cvu_hetero-'));
    growth=out.f;

    J = sum((G_data-growth).^2 );
 

