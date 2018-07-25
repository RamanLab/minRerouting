function [flux1,flux2,fluxWT] = myMOMA(model,rxn1,rxn2)
    
    rxn1=find(strcmp(rxn1,model.rxns));
    rxn2=find(strcmp(rxn2,model.rxns));
    
    solWT=optimizeCbModel(model,'max','zero'); 
    fluxWT=solWT.x;
    
    modeldel_1=model;
    modeldel_1.lb(rxn1)=0;
    modeldel_1.ub(rxn1)=0;
    
    modeldel_2=model;
    modeldel_2.lb(rxn2)=0;
    modeldel_2.ub(rxn2)=0;
    
    sol_11=fluxMOMAzn(modeldel_1,fluxWT); 
    flux1=sol_11.x;
    
    sol_21=fluxMOMAzn(modeldel_2,fluxWT);   
    flux2=sol_21.x;
    
end