function [flux1,flux2] = myMOMA(model,rxn1,rxn2)
    modeldel_1=model;
    modeldel_1.lb(rxn1)=0;
    modeldel_1.ub(rxn1)=0;
    
    modeldel_2=model;
    modeldel_2.lb(rxn2)=0;
    modeldel_2.ub(rxn2)=0;
    
    V1=zeros(length(model.rxns),2);
    V2=zeros(length(model.rxns),2);

    sol_11=optimizeCbModel(modeldel_1,'max','zero'); %Sol_11 is MT1 Iteration 1
    V1(:,1)=sol_11.x;
    
    sol_21=fluxMOMAzn(modeldel_2,V1(:,1));   %Sol_21 is MT2 Iteration 1
    V2(:,1)=sol_21.x;
    
    for j=2
        sol_1j=fluxMOMAzn(modeldel_1,V2(:,j-1));
        V1(:,j)=sol_1j.x;
        
        sol_2j=fluxMOMAzn(modeldel_2,V1(:,j));        
        V2(:,j)=sol_2j.x;
    end
    
    flux1=V1(:,2);
    flux2=V2(:,2);
    
end