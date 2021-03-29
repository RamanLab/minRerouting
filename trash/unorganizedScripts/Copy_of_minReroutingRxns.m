function [minRerouting] = minReroutingRxns(model,Jdl,cutOff,Iter)
%% [minRerouting] = minReroutingRxns(model,Jdl,cutOff)
% INPUT
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         Reaction Names
% Jdl            List of reaction pairs(double lethals generally) for identifying the flux rerouting
%OPTIONAL
% cutoff         cutoff flux difference value for MOMA difference.Default is 0.0001.
%
%OUTPUT
% minRerouting   The structure with reaction sets in alternate routes
%   rxns         List of total reactions in minimal reoruting set for each pair
%   diff         The flux difference value obtained after L0-MOMA
%   PathShort    List of reactions in shorter of the alternate paths
%   PathLong     List of reactions in longer of the alternate paths
%   PathCommon   List of reactions common in both the alternate paths
%
% Omkar Mohite       13 Jul,2017.
if (nargin <3 || isempty(cutOff))
        cutOff = 0.0001;
end

if (nargin <4 || isempty(Iter))
        Iter = 5;
end

[nLethals,temp] = size(Jdl);

% minRerouting consists information on reactions in alternate paths
minRerouting(nLethals).rxns=[];
minRerouting(nLethals).diff=[];
minRerouting(nLethals).flux1=[];
minRerouting(nLethals).flux2=[];

%%
h = waitbar(0,'0.00','Name','Identifying minRerouitngSets...');
modeldel_1=model;
modeldel_2=model;
for iLeth=1:nLethals
    delIdx_1=find(strcmp(Jdl(iLeth,1),model.rxns));
    delIdx_2=find(strcmp(Jdl(iLeth,2),model.rxns));
    
    modeldel_1.lb(delIdx_1)=0;
    modeldel_1.ub(delIdx_1)=0;
    
    modeldel_2.lb(delIdx_2)=0;
    modeldel_2.ub(delIdx_2)=0;
    
%     sol_11=optimizeCbModel(modeldel_1,'max','zero'); %Sol_11 is for modelDel1 -FBA
%     fluxFBA=sol_11.x;
%     
%     sol_21=fluxMOMAzn(modeldel_2,fluxFBA);   %Sol_21 is for modelDel2 close to V11 - MOMA
%     flux2=sol_21.x;
%     
%     sol_12=fluxMOMAzn(modeldel_1,flux2);   %Sol_12 is for modelDel1 close to V21 - MOMA
%     flux1=sol_12.x;
% %     
%     sol_22=fluxMOMAzn(modeldel_2,flux1);   %Sol_21 is for modelDel2 close to V11 - MOMA
%     flux2=sol_22.x;

% Run 3 iterations of MOMA
    V1=zeros(length(model.rxns),Iter);
    V2=zeros(length(model.rxns),Iter);

    sol_11=optimizeCbModel(modeldel_1,'max','zero'); %Sol_11 is MT1 Iteration 1
    V1(:,1)=sol_11.x;
    
    sol_21=fluxMOMAzn(modeldel_2,V1(:,1));   %Sol_21 is MT2 Iteration 1
    V2(:,1)=sol_21.x;
    
    for j=2:Iter
        sol_1j=fluxMOMAzn(modeldel_1,V2(:,j-1));
        V1(:,j)=sol_1j.x;
        
        sol_2j=fluxMOMAzn(modeldel_2,V1(:,j));        
        V2(:,j)=sol_2j.x;
    end
 
%  % Quadratic MOMA
%  sol_11=optimizeCbModel(modeldel_1,'max','one'); %Sol_11 is MT1 Iteration 1
%     V1(:,1)=sol_11.x;
%     
%     sol_21=fluxMOMA(modeldel_2,V1(:,1));   %Sol_21 is MT2 Iteration 1
%     V2(:,1)=sol_21.x;
%     
%     for j=2:3
%         sol_1j=fluxMOMA(modeldel_1,V2(:,j-1));
%         V1(:,j)=sol_1j.x;
%         
%         sol_2j=fluxMOMA(modeldel_2,V1(:,j));        
%         V2(:,j)=sol_2j.x;
%     end
    
    flux1=V1(:,Iter);
    flux2=V2(:,Iter);
    
    diff=flux1-flux2;
    
    minRerouting(iLeth).rxns=model.rxns(find(diff));
    minRerouting(iLeth).diff=diff(find(diff));
    minRerouting(iLeth).flux1=V1;
    minRerouting(iLeth).flux2=V2;
%     
%     flux1Rxn=model.rxns(find(flux1));
%     flux2Rxn=model.rxns(find(flux2));
%     Path2=minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns,flux1Rxn));
%     Path1=minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns,flux2Rxn));
%     
%     minRerouting(iLeth).pathCommon=Path1(ismember(Path1,Path2));
%     path1_Ex=Path1(~ismember(Path1,flux1Rxn));
%     path2_Ex=Path2(~ismember(Path2,flux2Rxn));
%     
%     if length(path1_Ex)<=length(path2_Ex)
%         minRerouting(iLeth).PathShort=path1_Ex;
%         minRerouting(iLeth).PathLong=path2_Ex;
%     else
%         minRerouting(iLeth).PathShort=path2_Ex;
%         minRerouting(iLeth).PathLong=path1_Ex;
%     end
%     
    modeldel_1.lb(delIdx_1)=model.lb(delIdx_1);
    modeldel_1.ub(delIdx_1)=model.ub(delIdx_1);
    
    modeldel_2.lb(delIdx_2)=model.lb(delIdx_2);
    modeldel_2.ub(delIdx_2)=model.ub(delIdx_2);
    waitbar(iLeth/nLethals,h,[num2str(round(iLeth*100/nLethals)) '% completed...']);
end

close(h);
end
