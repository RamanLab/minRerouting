function [ minReFBA ] = reroutingFBA( model, Jdl )
%REROUTINGFBA Summary of this function goes here
%   Detailed explanation goes here
nLethals=length(Jdl);
minReFBA(nLethals).rxns=[];
minReFBA(nLethals).cluster=[];
modeldel_1=model;
modeldel_2=model;
h = waitbar(0,'0.00','Name','Identifying minRerouitngSets...');

for iLeth=1:nLethals
    delIdx_1=find(strcmp(Jdl(iLeth,1),model.rxns));
    delIdx_2=find(strcmp(Jdl(iLeth,2),model.rxns));
    
    modeldel_1.lb(delIdx_1)=0;
    modeldel_1.ub(delIdx_1)=0;
    
    modeldel_2.lb(delIdx_2)=0;
    modeldel_2.ub(delIdx_2)=0;
    
    sol_1=optimizeCbModel(modeldel_1,'max','one'); %Sol_11 is MT1 Iteration 1
    flux1=sol_1.x;
    
    sol_2=optimizeCbModel(modeldel_2,'max','one'); %Sol_11 is MT1 Iteration 1
    flux2=sol_2.x;
    
    diff=flux1-flux2;
    
    minReFBA(iLeth).rxns=model.rxns(find(diff>0.0001));
    
    flux1Rxn=model.rxns(find(flux1>0.0001));
    flux2Rxn=model.rxns(find(flux2>0.0001));
    Path2=flux2Rxn(~ismember(flux2Rxn,flux1Rxn));
    Path1=flux1Rxn(~ismember(flux1Rxn,flux2Rxn));
    
    minReFBA(iLeth).cluster=[Path1;Path2];
    modeldel_1.lb(delIdx_1)=model.lb(delIdx_1);
    modeldel_1.ub(delIdx_1)=model.ub(delIdx_1);
    
    modeldel_2.lb(delIdx_2)=model.lb(delIdx_2);
    modeldel_2.ub(delIdx_2)=model.ub(delIdx_2);
    waitbar(iLeth/nLethals,h,[num2str(round(iLeth*100/nLethals)) '% completed...']);
end
close(h);
end