function [coreRerouting] = coreRerouting(model,Jdl)
[nLethals,temp] = size(Jdl);
cutOff = 0.000001;
% minRerouting consists information on reactions in alternate paths
coreRerouting(nLethals).rxns=[];
coreRerouting(nLethals).diff=[];
coreRerouting(nLethals).V1=[];
coreRerouting(nLethals).V2=[];

%%
h = waitbar(0,'0.00','Name','Identifying core-minRerouitngSets...');
modeldel_1=model;
modeldel_2=model;
for iLeth=1:nLethals
    delIdx_1=find(strcmp(Jdl(iLeth,1),model.rxns));
    delIdx_2=find(strcmp(Jdl(iLeth,2),model.rxns));
    
    modeldel_1.lb(delIdx_1)=0;
    modeldel_1.ub(delIdx_1)=0;
    
    modeldel_2.lb(delIdx_2)=0;
    modeldel_2.ub(delIdx_2)=0;
%% Run 3 iterations of MOMA
Iter=3;
    V1=zeros(length(model.rxns),Iter);
    V2=zeros(length(model.rxns),Iter);

    sol_11=optimizeCbModel(modeldel_1,'max','zero'); %Sol_11 is MT1 Iteration 1
    sol_11.x(sol_11.x<cutOff)=0;
    V1(:,1)=sol_11.x;
    
    sol_21=fluxMOMAzn(modeldel_2,V1(:,1));   %Sol_21 is MT2 Iteration 1
    sol_21.x(sol_21.x<cutOff)=0;
    V2(:,1)=sol_21.x;
    
    for j=2:Iter
        sol_1j=fluxMOMAzn(modeldel_1,V2(:,j-1));
        sol_1j.x(sol_1j.x<cutOff)=0;
        V1(:,j)=sol_1j.x;
        
        sol_2j=fluxMOMAzn(modeldel_2,V1(:,j));        
        sol_2j.x(sol_2j.x<cutOff)=0;
        V2(:,j)=sol_2j.x;
    end
    
    coreRerouting(nLethals).V1=V1;
    coreRerouting(nLethals).V2=V2;
    flux1=V1(:,Iter);
    flux2=V2(:,Iter);
    
    diff=flux1-flux2;
    min_ids=find((flux1==0 & flux2~=0) | (flux1~=0 & flux2==0));
    coreRerouting(iLeth).rxns=model.rxns(min_ids);
    coreRerouting(iLeth).diff=diff(min_ids);
    
    modeldel_1.lb(delIdx_1)=model.lb(delIdx_1);
    modeldel_1.ub(delIdx_1)=model.ub(delIdx_1);
    
    modeldel_2.lb(delIdx_2)=model.lb(delIdx_2);
    modeldel_2.ub(delIdx_2)=model.ub(delIdx_2);
    waitbar(iLeth/nLethals,h,[num2str(round(iLeth*100/nLethals)) '% completed...']);
end

close(h);
end
