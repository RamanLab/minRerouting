function [minRerouting] = minReroutingRxns(model,Jdl,cutOff)
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
if exist('cutoff', 'var')
    if isempty(cutoff)
        cutoff = 0.0001;
    end
else
    cutoff = 0.0001;
end

nLethals = length(Jdl);

% minRerouting consists information on reactions in alternate paths
minRerouting(nLethals).rxns=[];
minRerouting(nLethals).diff=[];
minRerouting(nLethals).PathShort=[];
minRerouting(nLethals).PathLong=[];
minRerouting(nLethals).pathCommon=[];

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
    
    sol_11=optimizeCbModel(modeldel_1,'max','zero'); %Sol_11 is for modelDel1 -FBA
    fluxFBA=sol_11.x;
    
    sol_21=fluxMOMAzn(modeldel_2,fluxFBA);   %Sol_21 is for modelDel2 close to V11 - MOMA
    flux2=sol_21.x;
    
    sol_12=fluxMOMAzn(modeldel_1,flux2);   %Sol_12 is for modelDel1 close to V21 - MOMA
    flux1=sol_12.x;

    diff=flux1-flux2;
    
    minRerouting(iLeth).rxns=model.rxns(find(abs(diff)>cutOff));
    minRerouting(iLeth).diff=diff(find(abs(diff)>cutOff));
    
    flux1Rxn=model.rxns(find(flux1));
    flux2Rxn=model.rxns(find(flux2));
    Path2=minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns,flux1Rxn));
    Path1=minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns,flux2Rxn));
    
    minRerouting(iLeth).pathCommon=Path1(ismember(Path1,Path2));
    path1_Ex=Path1(~ismember(Path1,flux1Rxn));
    path2_Ex=Path2(~ismember(Path2,flux2Rxn));
    
    if length(path1_Ex)<=length(path2_Ex)
        minRerouting(iLeth).PathShort=path1_Ex;
        minRerouting(iLeth).PathLong=path2_Ex;
    else
        minRerouting(iLeth).PathShort=path2_Ex;
        minRerouting(iLeth).PathLong=path1_Ex;
    end
    
    modeldel_1.lb(delIdx_1)=model.lb(delIdx_1);
    modeldel_1.ub(delIdx_1)=model.lb(delIdx_1);
    
    modeldel_2.lb(delIdx_2)=model.lb(delIdx_2);
    modeldel_2.ub(delIdx_2)=model.lb(delIdx_2);
    waitbar(iLeth/nLethals,h,[num2str(round(iLeth*100/nLethals)) '% completed...']);
end

close(h);
end
