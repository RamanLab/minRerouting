function [minRerouting] = minReroutingRxns(model,Jdl,cutOff,delta,Iter, Division)
%% [minRerouting] = minReroutingRxns(model,Jdl,cutOff)
% INPUT
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         Reaction Names
% Jdl            List of reaction pairs(doJdl(uble lethals generally) for identifying the flux rerouting
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
        cutOff = 0.000001; 
end

if (nargin <4 || isempty(delta))
        delta = 0.1; 
end

if (nargin <5 || isempty(Iter))
        Iter = 3;
end

if (nargin <6 || isempty(Division))
        Division = 'False';
end

[nLethals,temp] = size(Jdl);

% minRerouting consists information on reactions in alternate paths
minRerouting(nLethals).rxns=[];
minRerouting(nLethals).diff=[];
if strcmp(Division, 'True')
    minRerouting(nLethals).PathShort=[];
    minRerouting(nLethals).PathLong=[];
    minRerouting(nLethals).pathCommon=[];
end
%%
h = waitbar(0,'0.00','Name','Identifying minRerouitngSets...');
modeldel_1=model;
modeldel_2=model;
try
for iLeth=1:nLethals
    delIdx_1=find(strcmp(Jdl(iLeth,1),model.rxns));
    delIdx_2=find(strcmp(Jdl(iLeth,2),model.rxns));
    
    modeldel_1.lb(delIdx_1)=0;
    modeldel_1.ub(delIdx_1)=0;
    
    modeldel_2.lb(delIdx_2)=0;
    modeldel_2.ub(delIdx_2)=0;

%% Run 3 iterations of MOMA
    V1=zeros(length(model.rxns),Iter);
    V2=zeros(length(model.rxns),Iter);

    sol_11=optimizeCbModel(modeldel_1,'max','one'); %Sol_11 is MT1 Iteration 1
    sol_11.x(sol_11.x<cutOff)=0;
    V1(:,1)=sol_11.x;
    
    V1_stat(1,1)=sol_11.stat;
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
    
    flux1=V1(:,Iter);
    flux2=V2(:,Iter);
    
    diff=abs(flux1-flux2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    min_ids = find(diff>delta*abs(flux1) & diff>delta*abs(flux2) & diff>cutOff);

    minRerouting(iLeth).rxns=model.rxns(min_ids);
    minRerouting(iLeth).diff=diff(min_ids);
    
    if strcmp(Division, 'True')
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
    end
    
%     %% Testing Inconsistency 1 of number of iterations 
%     Diff_mat = abs(V1-V2);
%     minRerouting(iLeth).Diff1 = find(Diff_mat(:,1)>delta*V1(:,1) & Diff_mat(:,1)>delta*V2(:,1) & Diff_mat(:,1)>cutOff);
%     minRerouting(iLeth).Diff2 = find(Diff_mat(:,2)>delta*V1(:,2) & Diff_mat(:,2)>delta*V2(:,2) & Diff_mat(:,2)>cutOff);
%     minRerouting(iLeth).Diff3 = find(Diff_mat(:,3)>delta*V1(:,3) & Diff_mat(:,3)>delta*V2(:,3) & Diff_mat(:,3)>cutOff);
%     minRerouting(iLeth).Diff4 = find(Diff_mat(:,4)>delta*V1(:,4) & Diff_mat(:,4)>delta*V2(:,4) & Diff_mat(:,4)>cutOff);
%     minRerouting(iLeth).Diff5 = find(Diff_mat(:,5)>delta*V1(:,5) & Diff_mat(:,5)>delta*V2(:,5) & Diff_mat(:,5)>cutOff);
    
%% Check for Quadratic or L0Norm MOMA
%  sol_11=optimizeCbModel(modeldel_1,'max','one'); %Sol_11 is MT1 Iteration 1
%     V1(:,1)=sol_11.x;
%     sol_21=fluxMOMA(modeldel_2,V1(:,1));   %Sol_21 is MT2 Iteration 1
%     V2(:,1)=sol_21.x;     
%     for j=2:3
%         sol_1j=fluxMOMA(modeldel_1,V2(:,j-1));
%         V1(:,j)=sol_1j.x;
%         
%         sol_2j=fluxMOMA(modeldel_2,V1(:,j));        
%         V2(:,j)=sol_2j.x;
%     end
%%
    
    modeldel_1.lb(delIdx_1)=model.lb(delIdx_1);
    modeldel_1.ub(delIdx_1)=model.ub(delIdx_1);
    
    modeldel_2.lb(delIdx_2)=model.lb(delIdx_2);
    modeldel_2.ub(delIdx_2)=model.ub(delIdx_2);
    waitbar(iLeth/nLethals,h,[num2str(round(iLeth*100/nLethals)) '% completed...']);
end
end

close(h);
end
