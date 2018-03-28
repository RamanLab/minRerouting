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
    
    V1=zeros(length(model.rxns),3);
    V2=zeros(length(model.rxns),3);

    sol_21=optimizeCbModel(modeldel_2,'max','zero'); %Sol_11 is MT1 Iteration 1
    V2(:,1)=sol_21.x;
    
    sol_11=fluxMOMAzn(modeldel_1,V2(:,1));   %Sol_21 is MT2 Iteration 1
    V1(:,1)=sol_11.x;
    
    for j=2:3
        sol_2j=fluxMOMAzn(modeldel_2,V1(:,j-1));
        V2(:,j)=sol_2j.x;
        
        sol_1j=fluxMOMAzn(modeldel_1,V2(:,j));        
        V1(:,j)=sol_1j.x;
    end
    
    flux1=V1(:,3);
    flux2=V2(:,3);
    
end