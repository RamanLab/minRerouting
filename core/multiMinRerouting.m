function [Castle] = multiMinRerouting(Castle, minNorm, obj_slack)
% [Castle] = multiMinRerouting(Castle, minNorm) 
% multiMinRerouting identifies minimal rerouting sets
% for synthetic lethals of multiple models in Castle structure 
% Input:
%   Castle :    The structure with models and synthetic lethals data
%   minNorm :   'two', 'one' or 'zero' identifies minReroutingSets 
%               using minimal L2, L1 or L0 norm formulation of LP
% OPTIONAL
%   obj_slack   Permissible slack on the WT growth rate in the MoMA
%               default value: 0.05
%
% Omkar Satyavan Mohite 06 Aug, 2018
% N Sowmya Manojna   11 Oct, 2021

if (nargin < 3 || isempty(obj_slack))
    obj_slack = 0.05;
end

nModels = length(Castle.data);
for j= 1:nModels
    switch minNorm
        case 'two'
            disp("Finding set of minimal rerouting reactions using L2 minNorm in " + Castle.data(j).model_name);
            % quadraticMOMA version
            Castle.data(j).minRe = minReroutingRxns(Castle.data(j).model, Castle.data(j).Jdl, 2, obj_slack);
            fname = "results/" + Castle.data(j).model_name + "/" + "Castle_two_norm.mat";
        case 'one'
            disp("Finding set of minimal rerouting reactions using L1 minNorm in " + Castle.data(j).model_name);
            % linearMOMA version
            Castle.data(j).minRe = minReroutingRxns(Castle.data(j).model, Castle.data(j).Jdl, 1, obj_slack);
            fname = "results/" + Castle.data(j).model_name + "/" + "Castle_one_norm.mat";

        case 'zero'
            disp("Finding set of minimal rerouting reactions using L0 minNorm (iterative formulation) in " + Castle.data(j).model_name);
            % sparseMOMA version
            Castle.data(j).minRe = minReroutingRxns(Castle.data(j).model, Castle.data(j).Jdl, 0, obj_slack);
            fname = "results/" + Castle.data(j).model_name + "/" + "Castle_zero_norm.mat";
    end

    data = Castle.data(j);
    try
        save(fname, 'data')
    catch
        mkdir(strcat('results/', Castle.data(j).model_name))
        save(fname, 'data')
    end
end    
end


function [minRerouting] = minReroutingRxns(model, Jdl, p_norm, obj_slack, cutOff, delta, Division)
%% [minRerouting] = minReroutingRxns(model, Jdl, obj_slack, cutOff, delta, Division)
% INPUT
%   model            (the following fields are required - others can be supplied)
%       S            Stoichiometric matrix
%       b            Right hand side = dx/dt
%       c            Objective coefficients
%       lb           Lower bounds
%       ub           Upper bounds
%       rxns         Reaction Names
%   Jdl              List of reaction pairs(double lethals generally)
%                    for identifying the flux rerouting
% OPTIONAL
%   p_norm           Norm of optimization to be carried out (Default: 1)
%   obj_slack        Permissible slack on the WT growth rate in the MoMA (Default: 0.05)
%   cutOff           cutoff flux difference value for MOMA difference (Default is 0.0001)
%
% OUTPUT
%   minRerouting        The structure with reaction sets in alternate routes
%       rxns            List of total reactions in minimal rerouting set for each pair
%       diff_flux       The flux difference value obtained after P-Norm MOMA
%       abs_diff_flux   The abs flux difference value obtained after P-Norm MOMA
%       PathShort       List of reactions in shorter of the alternate paths
%       PathLong        List of reactions in longer of the alternate paths
%       PathCommon      List of reactions common in both the alternate paths
%
% N Sowmya Manojna   11 Oct, 2021

    if (nargin < 3 || isempty(p_norm))
            p_norm = 1;
    end

    if (nargin < 4 || isempty(obj_slack))
            obj_slack = 0.05;
    end

    if (nargin < 5 || isempty(cutOff))
            cutOff = 1e-6;
    end

    if (nargin < 6 || isempty(delta))
            delta = 0.05;
    end

    if (nargin < 7 || isempty(Division))
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

    % Finding the WT solution for grWT (Addition)
    sol = optimizeCbModel(model, 'max', 'one');
    grWT = sol.f;

    %%
    h = waitbar(0,'0.00','Name','Identifying minRerouitngSets...');
    modelDel1 = model;
    modelDel2 = model;
    for iLeth = 1:nLethals
        delIdx_1 = find(strcmp(Jdl(iLeth,1), model.rxns));
        delIdx_2 = find(strcmp(Jdl(iLeth,2), model.rxns));

        modelDel1.lb(delIdx_1) = 0;
        modelDel1.ub(delIdx_1) = 0;

        modelDel2.lb(delIdx_2) = 0;
        modelDel2.ub(delIdx_2) = 0;

        % fprintf('Finding minimal rerouting for pair: ( %s , %s )', Jdl(iLeth,1),Jdl(iLeth,2));

        switch p_norm
            case 0
                % Run the Sparse MoMA for the two deletions.
                [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = sparseMOMA_doubleKO(modelDel1, modelDel2, obj_slack, 'max');
            case 1
                % Run the Linear MoMA for the two deletions.
                [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = linearMOMA_doubleKO(modelDel1, modelDel2, obj_slack, 'max');
            case 2
                % Run the Quadratic MoMA for the two deletions.
                [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = quadraticMOMA_doubleKO(modelDel1, modelDel2, obj_slack, 'max');
        end

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
            % default: delta = 0.05;
            min_ids = find(diff_flux>delta*abs(flux1) & diff_flux>delta*abs(flux2) & diff_flux>cutOff);

            minRerouting(iLeth).del_rxn1 = del_rxn1;
            minRerouting(iLeth).del_rxn2 = del_rxn2;
            minRerouting(iLeth).solStatus = solStatus;
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
            minRerouting(iLeth).diff_flux = [];
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


function [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = sparseMOMA_doubleKO(modelDel1, modelDel2, obj_slack, osenseStr, verbFlag)
% Performs a sparse version of the MOMA (minimization of metabolic
% adjustment) approach upgraded for comparing double knockouts
%
% USAGE:
%    [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = sparseMOMA_doubleKO(modelDel1, modelDel2, osenseStr, minFluxFlag, verbFlab)
%
% INPUTS:
%    modelDe11:         Deletion strain model for Rxn1
%    modelDel2:         Deletion strain model for Rxn2
%
% OPTIONAL INPUTS:
%    osenseStr:        Maximize ('max') / minimize ('min') (Default = 'max')
%    verbFlag:         Verbose output (Default = false)
%
% OUTPUTS:
%    solutionDel1:     Deletion 1 solution structure
%    solutionDel2:     Deletion 2 solution structure
%    totalFluxDiff:    Value of the linear MOMA objective, i.e. :math:`\sum |v_{del1}-v_{del2}|`
%    solStatus:        Solution status - solves the problem: (`f_wt` is the optimal wild type objective value found by FBA)
%
% .. math::
%     min ~&~  \sum |v_{del1} - v_{del2}| \\
%         ~&~ S_{del1}v_{del1} = 0 \\
%         ~&~ lb_{del1} \leq v_{del1} \leq ub_{del1} \\
%         ~&~ c_{del1}^T v_{del1} = f_{del1} \\
%         ~&~ S_{del2}v_{del2} = 0 \\
%         ~&~ lb_{del2} \leq v_{del2} \leq ub_{del2}
%         ~&~ c_{del2}^T v_{del2} = f_{del2} \\
%%

% NOTE:
%
%    1) This formulation allows for selecting the most appropriate
%    optimal FBA solutions for both knockout strains as the starting point as opposed to
%    picking an arbitrary starting point (original MOMA implementation).
%
%    2) The reaction sets in the two models do not have to be equal as long as
%    there is at least one reaction in common
%
% .. Author: - Markus Herrgard 11/7/06
% ..           Omkar Satyavan Mohite 06/08/2018
% ..           N Sowmya Manojna 11/10/2021

    if (nargin < 3 || isempty(obj_slack))
        obj_slack = 0.05;
    end

    if (nargin < 4 || isempty(osenseStr))
        osenseStr = 'max';
    end

    if (nargin < 5)
        verbFlag = false;
    end

    % LP solution tolerance
    global CBT_LP_PARAMS
    if (exist('CBT_LP_PARAMS', 'var'))
        if isfield(CBT_LP_PARAMS, 'objTol')
            tol = CBT_LP_PARAMS.objTol;
        else
            tol = 1e-6;
        end
    else
        tol = 1e-6;
    end

    [nMets1,nRxns1] = size(modelDel1.S);
    [nMets2,nRxns2] = size(modelDel2.S);

    % Match model reaction sets
    commonRxns = ismember(modelDel1.rxns, modelDel2.rxns);
    nCommon = sum(commonRxns);
    if (nCommon == 0)
        error('No common rxns in the models');
    end

    solutionDel1.f = [];
    solutionDel1.x = [];
    solutionDel1.stat = -1;
    solutionDel2.f = [];
    solutionDel2.x = [];
    solutionDel2.stat = -1;

    if (verbFlag)
        fprintf('Solving FBA for deletion 1 strain: %d constraints %d variables ', nMets1, nRxns1);
    end

    % Solve wt problem
    solutionDel1 = optimizeCbModel(modelDel1, osenseStr, 'one');

    if (verbFlag)
        fprintf('%f seconds\n', solutionDel1.time);
    end

    % Round off solution to avoid numerical problems
    if (strcmp(osenseStr,'max'))
        objValDel1 = floor(solutionDel1.f/tol)*tol;
    else
        objValDel1 = ceil(solutionDel1.f/tol)*tol;
    end

    if (verbFlag)
        fprintf('Solving FBA for deletion 2 strain: %d constraints %d variables ', nMets1, nRxns2);
    end

    % Solve wt problem
    solutionDel2 = optimizeCbModel(modelDel2, osenseStr, 'one');

    if (verbFlag)
        fprintf('%f seconds\n', solutionDel2.time);
    end

    % Round off solution to avoid numerical problems
    if (strcmp(osenseStr,'max'))
        objValDel2 = floor(solutionDel2.f/tol)*tol;
    else
        objValDel2 = ceil(solutionDel2.f/tol)*tol;
    end

    % Variables in the following problem are
    % x = [v1;v2;delta+;delta-]
    % where v1 = deletion strain 1 flux vector
    %       v2 = deletion strain 2 flux vector
    %       delta+ = v1 - v2
    %       delta- = v2 - v1
    % delta+ and delta- constraints help reducing the common fluxes

    if (solutionDel1.stat > 0 && solutionDel2.stat > 0)
        % Construct the LHS matrix
        % Rows:
        % 1: Sdel1*v1 = 0 for the deletion strain 1
        % 2: Sdel2*v2 = 0 for the deletion strain 2
        % 3: delta+ >= v1-v2
        % 4: delta- >= v2-v1
        % 5: c'v1 >= 0.9*f1 (deletion strain 1) (10% slack on obj)
        % 6: c'v2 >= 0.9*f2 (deletion strain 2)
        % OR 5,6 : c'v1 and c'v2 >= 0.05*grWT
        % obj_slack = 0.1;

        A = [modelDel1.S sparse(nMets1,nRxns2+2*nCommon);
             sparse(nMets2,nRxns1) modelDel2.S sparse(nMets2,2*nCommon);
             createDeltaMatchMatrix(modelDel1.rxns,modelDel2.rxns);
             modelDel1.c' sparse(1,nRxns2+2*nCommon);
             sparse(1,nRxns1) modelDel2.c' sparse(1,2*nCommon)];

        % Construct the RHS vector
        b = [zeros(nMets1+nMets2+2*nCommon,1); (1-obj_slack)*objValDel1; (1-obj_slack)*objValDel2];

        % Construct the objective (sum of all delta+ and delta-)
        c = [zeros(nRxns1+nRxns2,1); ones(2*nCommon,1)];

        % Construct the ub/lb
        % delta+ and delta- are in [0 10000]
        lb = [modelDel1.lb; modelDel2.lb; zeros(2*nCommon,1)];
        ub = [modelDel1.ub; modelDel2.ub; 10000*ones(2*nCommon,1)];

        % Construct the constraint direction vector (G for delta's, E for
        % everything else)
        csense(1:(nMets1+nMets2)) = 'E';
        csense((nMets1+nMets2)+1:(nMets1+nMets2+2*nCommon)) = 'G';
        if (strcmp(osenseStr,'max'))
            csense(end+1) = 'G';
            csense(end+1) = 'G';
        else
            csense(end+1) = 'L';
            csense(end+1) = 'L';
        end

        if (verbFlag)
            fprintf('Solving 0-norm MOMA for double knockouts: %d constraints %d variables ',size(A,1),size(A,2));
        end

        [LPproblem.A, LPproblem.b, LPproblem.c, LPproblem.lb, LPproblem.ub, LPproblem.csense, LPproblem.osense] = deal(A, b, c, lb, ub, csense, 1);

        % Added cplex_direct support for zero norm
        % In order to use cplex_direct without any problems,
        % Consider cloning the latest version of cobratoolbox
        % The past versions have a problem with the solver initialization.
        % The autors of this paper raised a PR (now merged) that solves the problem.
        changeCobraSolver('cplex_direct', 'LP');
        LPsolution = solveCobraLP(LPproblem, 'minNorm', 0);
        changeCobraSolver('gurobi', 'LP');

        if (verbFlag)
            fprintf('%f seconds\n', LPsolution.time);
        end

        if (LPsolution.stat > 0)
            solutionDel1.x = LPsolution.full(1:nRxns1);
            solutionDel1.f = sum(modelDel1.c.*solutionDel1.x);
            solutionDel2.x = LPsolution.full((nRxns1+1):(nRxns1+nRxns2));
            solutionDel2.f = sum(modelDel2.c.*solutionDel2.x);
            totalFluxDiff = LPsolution.obj;
        else
            warning('sparse MOMA problem is infeasible or unconstrained');
            solStatus = LPsolution.stat;
            totalFluxDiff = [];
        end

        solutionDel1.stat = LPsolution.stat;
        solutionDel2.stat = LPsolution.stat;
        solStatus = LPsolution.stat;

    elseif solutionDel1.stat <= 0
        warning('Deletion 1 strain FBA problem is infeasible or unconstrained');
        solStatus = solutionDel1.stat;
        totalFluxDiff = [];
    elseif solutionDel2.stat <= 0
        totalFluxDiff = [];
        warning('Deletion 2 strain FBA problem is infeasible or unconstrained');
        solStatus = solutionDel2.stat;
    end
end

function [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = linearMOMA_doubleKO(modelDel1, modelDel2, obj_slack, osenseStr, verbFlag)
% Performs a linear version of the MOMA (minimization of metabolic
% adjustment) approach upgraded for comparing double knockouts
%
% USAGE:
%
%    [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = linearMOMA_doubleKO(modelDel1, modelDel2, osenseStr, minFluxFlag, verbFlab)
%
% INPUTS:
%    modelDe11:        Deletion strain model for Rxn1
%    modelDel2:        Deletion strain model for Rxn2
%
% OPTIONAL INPUTS:
%    obj_slack         Permissible slack on the WT growth rate in the MoMA (Default: 0.05)
%    osenseStr:        Maximize ('max') / minimize ('min') (Default = 'max')
%    verbFlag:         Verbose output (Default = false)
%
% OUTPUTS:
%    solutionDel1:     Deletion 1 solution structure
%    solutionDel2:     Deletion 2 solution structure
%    totalFluxDiff:    Value of the linear MOMA objective, i.e. :math:`\sum |v_{del1}-v_{del2}|`
%    solStatus:        Solution status - solves the problem: (`f_wt` is the optimal wild type objective value found by FBA)
%
% .. math::
%     min ~&~  \sum |v_{del1} - v_{del2}| \\
%         ~&~ S_{del1}v_{del1} = 0 \\
%         ~&~ lb_{del1} \leq v_{del1} \leq ub_{del1} \\
%         ~&~ c_{del1}^T v_{del1} = f_{del1} \\
%         ~&~ S_{del2}v_{del2} = 0 \\
%         ~&~ lb_{del2} \leq v_{del2} \leq ub_{del2}
%         ~&~ c_{del2}^T v_{del2} = f_{del2} \\
%
% NOTE:
%
%    1) This formulation allows for selecting the most appropriate
%    optimal FBA solutions for both knockout strains as the starting point as opposed to
%    picking an arbitrary starting point (original MOMA implementation).
%
%    2) The reaction sets in the two models do not have to be equal as long as
%    there is at least one reaction in common
%
% .. Author: - Markus Herrgard 11/7/06
% ..           Omkar Satyavan Mohite 06/08/2018

    if (nargin < 3 || isempty(obj_slack))
        obj_slack = 0.05;
    end

    if (nargin < 4 || isempty(osenseStr))
        osenseStr = 'max';
    end

    if (nargin < 5 || isempty(minFluxFlag))
        minFluxFlag = false;
    end

    if (nargin < 6)
        verbFlag = false;
    end

    % LP solution tolerance
    global CBT_LP_PARAMS
    if (exist('CBT_LP_PARAMS', 'var'))
        if isfield(CBT_LP_PARAMS, 'objTol')
            tol = CBT_LP_PARAMS.objTol;
        else
            tol = 1e-6;
        end
    else
        tol = 1e-6;
    end

    [nMets1,nRxns1] = size(modelDel1.S);
    [nMets2,nRxns2] = size(modelDel2.S);

    % Match model reaction sets
    commonRxns = ismember(modelDel1.rxns, modelDel2.rxns);
    nCommon = sum(commonRxns);
    if (nCommon == 0)
        error('No common rxns in the models');
    end

    solutionDel1.f = [];
    solutionDel1.x = [];
    solutionDel1.stat = -1;
    solutionDel2.f = [];
    solutionDel2.x = [];
    solutionDel2.stat = -1;

    if (verbFlag)
        fprintf('Solving FBA for deletion 1 strain: %d constraints %d variables ', nMets1, nRxns1);
    end
    % Solve wt problem
    solutionDel1 = optimizeCbModel(modelDel1, osenseStr, 'one');

    if (verbFlag)
        fprintf('%f seconds\n', solutionDel1.time);
    end
    % Round off solution to avoid numerical problems
    if (strcmp(osenseStr,'max'))
        objValDel1 = floor(solutionDel1.f/tol)*tol;
    else
        objValDel1 = ceil(solutionDel1.f/tol)*tol;
    end

    if (verbFlag)
        fprintf('Solving FBA for deletion 2 strain: %d constraints %d variables ', nMets1, nRxns2);
    end
    % Solve wt problem
    solutionDel2 = optimizeCbModel(modelDel2, osenseStr, 'one');

    if (verbFlag)
        fprintf('%f seconds\n', solutionDel2.time);
    end
    % Round off solution to avoid numerical problems
    if (strcmp(osenseStr,'max'))
        objValDel2 = floor(solutionDel2.f/tol)*tol;
    else
        objValDel2 = ceil(solutionDel2.f/tol)*tol;
    end

    % Variables in the following problem are
    % x = [v1;v2;delta+;delta-]
    % where v1 = deletion strain 1 flux vector
    %       v2 = deletion strain 2 flux vector
    %       delta+ = v1 - v2
    %       delta- = v2 - v1

    if (solutionDel1.stat > 0 && solutionDel2.stat > 0)
        % Construct the LHS matrix
        % Rows:
        % 1: Sdel1*v1 = 0 for the deletion strain 1
        % 2: Sdel2*v2 = 0 for the deletion strain 2
        % 3: delta+ >= v1-v2
        % 4: delta- >= v2-v1
        % 5: c'v1 >= 0.9*f1 (deletion strain 1) (10 % slack on obj)
        % 6: c'v2 >= 0.9*f2 (deletion strain 2)
        % OR 5,6: c'v1 and c'v2 >= 0.05*grWT
        % Conditions 3 and 4 ensure that the flux through
        % the common reactions are minimized.

        A = [modelDel1.S sparse(nMets1,nRxns2+2*nCommon);
             sparse(nMets2,nRxns1) modelDel2.S sparse(nMets2,2*nCommon);
             createDeltaMatchMatrix(modelDel1.rxns,modelDel2.rxns);
             modelDel1.c' sparse(1,nRxns2+2*nCommon);
             sparse(1,nRxns1) modelDel2.c' sparse(1,2*nCommon)];

        % Construct the RHS vector
        b = [zeros(nMets1+nMets2+2*nCommon,1); (1-obj_slack)*objValDel1; (1-obj_slack)*objValDel2];

        % Construct the objective (sum of all delta+ and delta-)
        c = [zeros(nRxns1+nRxns2,1); ones(2*nCommon,1)];

        % Construct the ub/lb
        % delta+ and delta- are in [0 10000]
        lb = [modelDel1.lb; modelDel2.lb; zeros(2*nCommon,1)];
        ub = [modelDel1.ub; modelDel2.ub; 10000*ones(2*nCommon,1)];

        % Construct the constraint direction vector (G for delta's, E for
        % everything else)
        csense(1:(nMets1+nMets2)) = 'E';
        csense((nMets1+nMets2)+1:(nMets1+nMets2+2*nCommon)) = 'G';
        if (strcmp(osenseStr,'max'))
            csense(end+1) = 'G';
            csense(end+1) = 'G';
        else
            csense(end+1) = 'L';
            csense(end+1) = 'L';
        end

        if (verbFlag)
            fprintf('Solving linear MOMA upgraded for double knockouts: %d constraints %d variables ',size(A,1),size(A,2));
        end

        % Solve the linearMOMA problem
        [LPproblem.A, LPproblem.b, LPproblem.c, LPproblem.lb, LPproblem.ub, LPproblem.csense, LPproblem.osense] = deal(A, b, c, lb, ub, csense, 1);
        LPsolution = solveCobraLP(LPproblem);

        if (verbFlag)
            fprintf('%f seconds\n', LPsolution.time);
        end

        if (LPsolution.stat > 0)
            solutionDel1.x = LPsolution.full(1:nRxns1);
            solutionDel1.f = sum(modelDel1.c.*solutionDel1.x);
            solutionDel2.x = LPsolution.full((nRxns1+1):(nRxns1+nRxns2));
            solutionDel2.f = sum(modelDel2.c.*solutionDel2.x);
            totalFluxDiff = LPsolution.obj;
        else
            warning('linear MOMA problem is infeasible or unconstrained');
            solStatus = LPsolution.stat;
            totalFluxDiff = [];
        end

        solutionDel1.stat = LPsolution.stat;
        solutionDel2.stat = LPsolution.stat;
        solStatus = LPsolution.stat;

    elseif solutionDel1.stat <= 0
        warning('Deletion 1 strain FBA problem is infeasible or unconstrained');
        solStatus = solutionDel1.stat;
        totalFluxDiff = [];
    elseif solutionDel2.stat <= 0
        totalFluxDiff = [];
        warning('Deletion 2 strain FBA problem is infeasible or unconstrained');
        solStatus = solutionDel2.stat;
    end
end

function [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = quadraticMOMA_doubleKO(modelDel1, modelDel2, obj_slack, osenseStr, verbFlag)
% Performs a quadratic version of the MOMA (minimization of metabolic
% adjustment) approach upgraded for comparing double knockouts
%
% USAGE:
%    [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = quadraticMOMA_doubleKO(modelDel1, modelDel2, osenseStr, minFluxFlag, verbFlab)
%
% INPUTS:
%    modelDe11:         Deletion strain model for Rxn1
%    modelDel2:         Deletion strain model for Rxn2
%
% OPTIONAL INPUTS:
%    osenseStr:        Maximize ('max') / minimize ('min') (Default = 'max')
%    verbFlag:         Verbose output (Default = false)
%
% OUTPUTS:
%    solutionDel1:     Deletion 1 solution structure
%    solutionDel2:     Deletion 2 solution structure
%    totalFluxDiff:    Value of the linear MOMA objective, i.e. :math:`\sum |v_{del1}-v_{del2}|`
%    solStatus:        Solution status - solves the problem: (`f_wt` is the optimal wild type objective value found by FBA)
%
% .. math::
%     min ~&~  \sum |v_{del1} - v_{del2}|_2 \\
%         ~&~ S_{del1}v_{del1} = 0 \\
%         ~&~ lb_{del1} \leq v_{del1} \leq ub_{del1} \\
%         ~&~ c_{del1}^T v_{del1} = f_{del1} \\
%         ~&~ S_{del2}v_{del2} = 0 \\
%         ~&~ lb_{del2} \leq v_{del2} \leq ub_{del2}
%         ~&~ c_{del2}^T v_{del2} = f_{del2} \\
%
% NOTE:
%
%    1) This formulation allows for selecting the most appropriate
%    optimal FBA solutions for both knockout strains as the starting point as opposed to
%    picking an arbitrary starting point (original MOMA implementation).
%
%    2) The reaction sets in the two models do not have to be equal as long as
%    there is at least one reaction in common
%
% .. Author: - Markus Herrgard 11/7/06
% ..           Omkar Satyavan Mohite 06/08/2018
% ..           N Sowmya Manojna 11/10/2021

    if (nargin < 3 || isempty(obj_slack))
        obj_slack = 0.05;
    end

    if (nargin < 4 || isempty(osenseStr))
        osenseStr = 'max';
    end

    if (nargin < 5)
        verbFlag = false;
    end

    % LP solution tolerance
    global CBT_LP_PARAMS
    if (exist('CBT_LP_PARAMS', 'var'))
        if isfield(CBT_LP_PARAMS, 'objTol')
            tol = CBT_LP_PARAMS.objTol;
        else
            tol = 1e-6;
        end
    else
        tol = 1e-6;
    end

    [nMets1,nRxns1] = size(modelDel1.S);
    [nMets2,nRxns2] = size(modelDel2.S);

    % Match model reaction sets
    commonRxns = ismember(modelDel1.rxns,modelDel2.rxns);
    nCommon = sum(commonRxns);
    if (nCommon == 0)
        error('No common rxns in the models');
    end

    solutionDel1.f = [];
    solutionDel1.x = [];
    solutionDel1.stat = -1;
    solutionDel2.f = [];
    solutionDel2.x = [];
    solutionDel2.stat = -1;

    if (verbFlag)
        fprintf('Solving FBA for deletion 1 strain: %d constraints %d variables ', nMets1, nRxns1);
    end

    % Solve wt problem
    solutionDel1 = optimizeCbModel(modelDel1,osenseStr, 'one');

    if (verbFlag)
        fprintf('%f seconds\n',solutionDel1.time);
    end

    % Round off solution to avoid numerical problems
    if (strcmp(osenseStr,'max'))
        objValDel1 = floor(solutionDel1.f/tol)*tol;
    else
        objValDel1 = ceil(solutionDel1.f/tol)*tol;
    end

    if (verbFlag)
        fprintf('Solving FBA for deletion 2 strain: %d constraints %d variables ',nMets1,nRxns2);
    end

    % Solve wt problem
    solutionDel2 = optimizeCbModel(modelDel2, osenseStr, 'one');

    if (verbFlag)
        fprintf('%f seconds\n', solutionDel2.time);
    end

    % Round off solution to avoid numerical problems
    if (strcmp(osenseStr, 'max'))
        objValDel2 = floor(solutionDel2.f/tol)*tol;
    else
        objValDel2 = ceil(solutionDel2.f/tol)*tol;
    end

    % Variables in the following problem are
    % x = [v1;v2;delta]
    % where v1 = wild type flux vector
    %       v2 = deletion strain flux vector
    %       delta = v1 - v2

    if (solutionDel1.stat > 0 && solutionDel2.stat > 0)
        % Construct the LHS matrix
        % Rows:
        % 1: Sdel1*v1 = 0 for the deletion strain 1
        % 2: Sdel2*v2 = 0 for the deletion strain 2
        % 3: delta+ >= v1-v2
        % 4: delta- >= v2-v1
        % 5: c'v1 >= 0.9*f1 (deletion strain 1) (10 % slack on obj)
        % 6: c'v2 >= 0.9*f2 (deletion strain 2)

        deltaMat = createDeltaMatchMatrix(modelDel1.rxns,modelDel2.rxns);
        deltaMat = deltaMat(1:nCommon,1:(nRxns1+nRxns2+nCommon));
        A = [modelDel1.S sparse(nMets1,nRxns2+nCommon);
            sparse(nMets2,nRxns1) modelDel2.S sparse(nMets2,nCommon);
            deltaMat;
            modelDel1.c' sparse(1,nRxns2+nCommon);
            sparse(1,nRxns1) modelDel2.c' sparse(1,nCommon)];

        % Construct the RHS vector
        b = [zeros(nMets1+nMets2+nCommon,1);(1-obj_slack)*objValDel1 ;(1-obj_slack)*objValDel2];

        c = [zeros(nRxns1+nRxns2+nCommon,1)];

        % Construct the ub/lb
        % delta [-10000 10000]
        lb = [modelDel1.lb;modelDel2.lb;-10000*ones(nCommon,1)];
        ub = [modelDel1.ub;modelDel2.ub;10000*ones(nCommon,1)];

        % Construct the constraint direction vector (G for delta's, E for
        % everything else)
        csense(1:(nMets1+nMets2+nCommon)) = 'E';
        if (strcmp(osenseStr,'max'))
            csense(end+1) = 'G';
            csense(end+1) = 'G';
        else
            csense(end+1) = 'L';
            csense(end+1) = 'L';
        end

        % F matrix
        F = [sparse(nRxns1+nRxns2,nRxns1+nRxns2+nCommon);
            sparse(nCommon,nRxns1+nRxns2) 2*eye(nCommon)];

        if (verbFlag)
            fprintf('Solving quadratic MOMA upgraded for double knockouts: %d constraints %d variables ',size(A,1),size(A,2));
        end

        % Solve the quadraticMOMA problem
        [QPproblem.A,QPproblem.b,QPproblem.F,QPproblem.c,QPproblem.lb,QPproblem.ub,QPproblem.csense,QPproblem.osense] = deal(A,b,F,c,lb,ub,csense,1);
        % QPsolution = solveCobraQP(QPproblem,[],verbFlag-1);
        QPsolution = solveCobraQP(QPproblem);

        if (verbFlag)
            fprintf('%f seconds\n',QPsolution.time);
        end

        % Get the solution(s)
        if QPsolution.stat == 1
            solutionDel1.x = QPsolution.full(1:nRxns1);
            solutionDel2.x = QPsolution.full((nRxns1+1):(nRxns1+nRxns2));

            solutionDel1.f = sum(modelDel2.c.*solutionDel1.x);
            solutionDel2.f = sum(modelDel2.c.*solutionDel2.x);

            totalFluxDiff = sum((solutionDel1.x-solutionDel2.x).^2);
        else
            totalFluxDiff = [];
        end
        solutionDel1.stat = QPsolution.stat;
        solutionDel2.stat = QPsolution.stat;
        solStatus = QPsolution.stat;

    elseif solutionDel1.stat <= 0
        warning('Deletion 1 strain FBA problem is infeasible or unconstrained');
        solStatus = solutionDel1.stat;
        totalFluxDiff = [];
    else
        warning('Deletion 2 strain FBA problem is infeasible or unconstrained');
        solStatus = solutionDel2.stat;
        totalFluxDiff = [];
    end
end
