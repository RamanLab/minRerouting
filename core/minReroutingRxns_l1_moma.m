function [minRerouting] = minReroutingRxns_l1_moma(model,Jdl,cutOff,delta,Division)
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
% OPTIONAL
% cutoff         cutoff flux difference value for MOMA difference.Default is 0.0001.
%
% OUTPUT
% minRerouting   The structure with reaction sets in alternate routes
%   rxns         List of total reactions in minimal reoruting set for each pair
%   diff         The flux difference value obtained after L0-MOMA
%   PathShort    List of reactions in shorter of the alternate paths
%   PathLong     List of reactions in longer of the alternate paths
%   PathCommon   List of reactions common in both the alternate paths
%
% Omkar Mohite       13 Jul,2017.

if (nargin < 3 || isempty(cutOff))
        cutOff = 0.000001; 
end

if (nargin < 4 || isempty(delta))
        delta = 0.05; 
end

if (nargin < 5 || isempty(Division))
        Division = 'True';
end

[nLethals,temp] = size(Jdl);

% minRerouting consists information on reactions in alternate paths
minRerouting(nLethals).rxns = [];
minRerouting(nLethals).diff = [];
if strcmp(Division, 'True')
    minRerouting(nLethals).PathShort = [];
    minRerouting(nLethals).PathLong = [];
    minRerouting(nLethals).pathCommon = [];
end
%%
solWT = optimizeCbModel(model,'max','one');
grWT = solWT.f;
h = waitbar(0,'0.00','Name','Identifying minRerouitngSets...');
modelDel1 = model;
modelDel2 = model;
for iLeth = 1:nLethals
    delIdx_1 = find(strcmp(Jdl(iLeth,1),model.rxns));
    delIdx_2 = find(strcmp(Jdl(iLeth,2),model.rxns));
    
    modelDel1.lb(delIdx_1) = 0;
    modelDel1.ub(delIdx_1) = 0;
    
    modelDel2.lb(delIdx_2) = 0;
    modelDel2.ub(delIdx_2) = 0;
    
    % fprintf('Finding minimal rerouting for pair: ( %s , %s )', Jdl(iLeth,1),Jdl(iLeth,2));
    
    [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = linearMOMA_doubleKO(modelDel1, modelDel2, grWT,'max');
    
    if solStatus > 0
        flux1 = solutionDel1.x;
        flux2 = solutionDel2.x;
        
        diff = abs(flux1-flux2);

        min_ids = find(diff>delta*abs(flux1) & diff>delta*abs(flux2) & diff>cutOff);

        minRerouting(iLeth).rxns = model.rxns(min_ids);
        minRerouting(iLeth).diff = diff(min_ids);
        minRerouting(iLeth).totalFluxDiff = totalFluxDiff;
    
        if strcmp(Division, 'True')
            flux1Rxn = model.rxns(find(flux1));
            flux2Rxn = model.rxns(find(flux2));
            Path2 = minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns, flux1Rxn));
            Path1 = minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns, flux2Rxn));

            minRerouting(iLeth).pathCommon = Path1(ismember(Path1, Path2));
            path1_Ex = Path1(~ismember(Path1, flux1Rxn));
            path2_Ex = Path2(~ismember(Path2, flux2Rxn));

            if length(path1_Ex) <= length(path2_Ex)
                minRerouting(iLeth).PathShort = path1_Ex;
                minRerouting(iLeth).PathLong = path2_Ex;
            else
                minRerouting(iLeth).PathShort = path2_Ex;
                minRerouting(iLeth).PathLong = path1_Ex;
            end
            
        end
    else
        minRerouting(iLeth).rxns = [];
        minRerouting(iLeth).diff = [];
        minRerouting(iLeth).totalFluxDiff = [];
        if strcmp(Division, 'True')
            minRerouting(iLeth).pathCommon = [];           
            minRerouting(iLeth).PathShort = [];
            minRerouting(iLeth).PathLong = [];
        end            
    end    
    
    minRerouting(iLeth).solStatus = solStatus;
    
    modelDel1.lb(delIdx_1) = model.lb(delIdx_1);
    modelDel1.ub(delIdx_1) = model.ub(delIdx_1);
    
    modelDel2.lb(delIdx_2) = model.lb(delIdx_2);
    modelDel2.ub(delIdx_2) = model.ub(delIdx_2);
    waitbar(iLeth/nLethals, h, [num2str(round(iLeth*100/nLethals)) '% completed...']);
end

close(h);
end
