function [minRerouting] = minReroutingRxns_l1_moma(model, Jdl, obj_slack, cutOff, delta, Division)
%% [minRerouting] = minReroutingRxns_l1_moma(model, Jdl, obj_slack, cutOff)
% INPUT
%   model (the following fields are required - others can be supplied)
%       S            Stoichiometric matrix
%       b            Right hand side = dx/dt
%       c            Objective coefficients
%       lb           Lower bounds
%       ub           Upper bounds
%       rxns         Reaction Names
%   Jdl              List of reaction pairs(double lethals generally)
%                    for identifying the flux rerouting
% OPTIONAL
%   obj_slack        Permissible slack on the WT growth rate in the MoMA (Default: 0.05)
%   cutOff           cutoff flux difference value for MOMA difference (Default is 0.0001)
%
% OUTPUT
%   minRerouting        The structure with reaction sets in alternate routes
%       rxns            List of total reactions in minimal rerouting set for each pair
%       diff_flux       The flux difference value obtained after L1-MOMA
%       abs_diff_flux   The abs flux difference value obtained after L1-MOMA
%       PathShort       List of reactions in shorter of the alternate paths
%       PathLong        List of reactions in longer of the alternate paths
%       PathCommon      List of reactions common in both the alternate paths
%
% Omkar Mohite       13 July, 2017
% N Sowmya Manojna   26 June, 2021

if (nargin < 3 || isempty(obj_slack))
        obj_slack = 0.05;
end

if (nargin < 4 || isempty(cutOff))
        cutOff = 1e-6; 
end

if (nargin < 5 || isempty(delta))
        delta = 0.05; 
end

if (nargin < 6 || isempty(Division))
        Division = 'True';
end

[nLethals, ~] = size(Jdl);

% minRerouting consists information on reactions in alternate paths
minRerouting(nLethals).rxns = [];
minRerouting(nLethals).diff_flux = [];
minRerouting(nLethals).abs_diff_flux = [];
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
    
    % Run the 1-norm Linear MoMA for the two deletions.
    [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = linearMOMA_doubleKO(modelDel1, modelDel2, obj_slack, 'max');
    
    if solStatus > 0
        flux1 = solutionDel1.x;
        flux2 = solutionDel2.x;

        % Track which two reactions were deleted
        del_rxn1 = Jdl(iLeth,1);
        del_rxn2 = Jdl(iLeth,2);
        
        
        % Get the flux difference between the two deletions
        diff_flux = flux1-flux2;
        abs_diff_flux = abs(flux1-flux2);

        % Get all the locations where the fluxes differ.
        % The two delta*abs() conditions ensure that the change in flux
        % is at least greater than a factor (delta) of the original flux.
        min_ids = find(abs_diff_flux>delta*abs(flux1) & abs_diff_flux>delta*abs(flux2) & abs_diff_flux>cutOff);

        % Find reactions that have modified fluxes
        minRerouting(iLeth).del_rxn1 = del_rxn1;
        minRerouting(iLeth).del_rxn2 = del_rxn2;
        minRerouting(iLeth).rxns = model.rxns(min_ids);
        minRerouting(iLeth).totalFluxDiff = totalFluxDiff;
        minRerouting(iLeth).diff_flux = diff_flux(min_ids);
        minRerouting(iLeth).abs_diff_flux = abs_diff_flux(min_ids);
    
        if strcmp(Division, 'True')
            % Get all active reactions in model 1 and 2
            flux1Rxn = model.rxns(find(flux1));
            flux2Rxn = model.rxns(find(flux2));

            % Find all reactions with different flux, that are active in model 1
            Path2 = minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns, flux1Rxn));
            % Find all reactions with different flux, that are active in model 2
            Path1 = minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns, flux2Rxn));
            % Find all reactions that are common between the two paths.
            minRerouting(iLeth).pathCommon = Path1(ismember(Path1, Path2));

            % Find exclusive differential reactions for each of the models
            path1_Ex = Path1(~ismember(Path1, flux1Rxn));
            path2_Ex = Path2(~ismember(Path2, flux2Rxn));

            % Assign long and short paths
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
        minRerouting(iLeth).abs_diff_flux = [];
        minRerouting(iLeth).totalFluxDiff = [];
        if strcmp(Division, 'True')
            minRerouting(iLeth).pathCommon = [];           
            minRerouting(iLeth).PathShort = [];
            minRerouting(iLeth).PathLong = [];
        end            
    end    
    
    minRerouting(iLeth).solStatus = solStatus;
    
    % Reset the LB and UB for next loop    
    modelDel1.lb(delIdx_1) = model.lb(delIdx_1);
    modelDel1.ub(delIdx_1) = model.ub(delIdx_1);
    
    modelDel2.lb(delIdx_2) = model.lb(delIdx_2);
    modelDel2.ub(delIdx_2) = model.ub(delIdx_2);
    waitbar(iLeth/nLethals, h, [num2str(round(iLeth*100/nLethals)) '% completed...']);
end

close(h);
end
