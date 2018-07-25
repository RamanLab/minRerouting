function [fluxDist] = doubleKO(model,Jdl)
%DOUBLEKO 
%   
modeldel_1=model;
modeldel_2=model;

for iLeth=1:length(Jdl)
    delIdx_1=find(strcmp(Jdl(iLeth,1),model.rxns));
    delIdx_2=find(strcmp(Jdl(iLeth,2),model.rxns));
    
    modeldel_1.lb(delIdx_1)=0;
    modeldel_1.ub(delIdx_1)=0;
    
    modeldel_2.lb(delIdx_2)=0;
    modeldel_2.ub(delIdx_2)=0;

    sol_KO1 = optimizeCbModel(modeldel_1,'max','zero'); %Sol_11 is MT1 Iteration 1
    sol_KO2 = optimizeCbModel(modeldel_2,'max','zero'); %Sol_11 is MT1 Iteration 1
    
    fluxDist(iLeth).KO1 = sol_KO1.x;
    fluxDist(iLeth).KO2 = sol_KO2.x;
end
end

