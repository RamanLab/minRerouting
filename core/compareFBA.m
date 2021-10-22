function [Castle] = compareFBA(Castle, minNorm, obj_slack)
    if (nargin < 3 || isempty(obj_slack))
        obj_slack = 0.05;
    end

    nModels = length(Castle.data);
    for j= 1:nModels
        switch minNorm
            case 'two'
                disp("Solving 2-norm FBA for the model " + Castle.data(j).model_name);
                % quadraticMOMA version
                Castle.data(j).fba_solution = fba_solution(Castle.data(j).model, Castle.data(j).Jdl, 2, obj_slack);
                fname = "results/" + Castle.data(j).model_name + "/" + "fba_two_norm.mat";
            case 'one'
                disp("Solving 1-norm FBA for the model " + Castle.data(j).model_name);
                % linearMOMA version
                Castle.data(j).fba_solution = fba_solution(Castle.data(j).model, Castle.data(j).Jdl, 1, obj_slack);
                fname = "results/" + Castle.data(j).model_name + "/" + "fba_one_norm.mat";

            case 'zero'
                disp("Solving 0-norm FBA for the model " + Castle.data(j).model_name);
                % sparseMOMA version
                Castle.data(j).fba_solution = fba_solution(Castle.data(j).model, Castle.data(j).Jdl, 0, obj_slack);
                fname = "results/" + Castle.data(j).model_name + "/" + "fba_zero_norm.mat";
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


function [serrano] = fba_solution(model, Jdl, p_norm, obj_slack, cutOff, delta, Division)
%% [serrano] = fba_solution(model, Jdl, obj_slack, cutOff, delta, Division)
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
%   serrano        The structure with reaction sets in alternate routes
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

    % serrano consists information on reactions in alternate paths
    serrano(nLethals).rxns = [];
    serrano(nLethals).diff_flux = [];
    serrano(nLethals).abs_diff_flux = [];

    if strcmp(Division, 'True')
        serrano(nLethals).PathShort = [];
        serrano(nLethals).PathLong = [];
        serrano(nLethals).pathCommon = [];
    end

    % Finding the WT solution for grWT (Addition)
    sol = optimizeCbModel(model, 'max', 'one');
    grWT = sol.f;

    %%
    h = waitbar(0,'0.00','Name','Identifying FBA Set...');
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
                [solutionDel1, solutionDel2, solStatus] = get_serrano_results(modelDel1, modelDel2, 0, obj_slack, 'max');
            case 1
                % Run the Linear MoMA for the two deletions.
                [solutionDel1, solutionDel2, solStatus] = get_serrano_results(modelDel1, modelDel2, 1, obj_slack, 'max');
            case 2
                % Run the Quadratic MoMA for the two deletions.
                [solutionDel1, solutionDel2, solStatus] = get_serrano_results(modelDel1, modelDel2, 2, obj_slack, 'max');
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

            serrano(iLeth).del_rxn1 = del_rxn1;
            serrano(iLeth).del_rxn2 = del_rxn2;
            serrano(iLeth).solStatus = solStatus;
            serrano(iLeth).rxns = model.rxns(min_ids);
            serrano(iLeth).diff_flux = diff_flux(min_ids);
            serrano(iLeth).abs_diff_flux = abs_diff_flux(min_ids);
            serrano(iLeth).totalFluxDiff = sum(abs_diff_flux(min_ids));

            if strcmp(Division, 'True')
                % Get all active reactions in model 1 and 2
                flux1Rxn = model.rxns(find(flux1));
                flux2Rxn = model.rxns(find(flux2));

                % Find all reactions with different flux, that are active in model 1
                Path2 = serrano(iLeth).rxns(ismember(serrano(iLeth).rxns, flux1Rxn));
                % Find all reactions with different flux, that are active in model 2
                Path1 = serrano(iLeth).rxns(ismember(serrano(iLeth).rxns, flux2Rxn));
                % Find all reactions that are common between the two paths.
                serrano(iLeth).pathCommon = Path1(ismember(Path1, Path2));

                % Find exclusive differential reactions for each of the models
                path1_Ex = Path1(~ismember(Path1, flux1Rxn));
                path2_Ex = Path2(~ismember(Path2, flux2Rxn));

                % Assign long and short paths
                if length(path1_Ex) <= length(path2_Ex)
                    serrano(iLeth).PathShort = path1_Ex;
                    serrano(iLeth).PathLong = path2_Ex;
                else
                    serrano(iLeth).PathShort = path2_Ex;
                    serrano(iLeth).PathLong = path1_Ex;
                end
            end
        else
            serrano(iLeth).rxns = [];
            serrano(iLeth).diff_flux = [];
            serrano(iLeth).abs_diff_flux = [];
            serrano(iLeth).totalFluxDiff = [];
            if strcmp(Division, 'True')
                serrano(iLeth).pathCommon = [];
                serrano(iLeth).PathShort = [];
                serrano(iLeth).PathLong = [];
            end
        end

        serrano(iLeth).solStatus = solStatus;

        % Reset the LB and UB for next loop
        modelDel1.lb(delIdx_1) = model.lb(delIdx_1);
        modelDel1.ub(delIdx_1) = model.ub(delIdx_1);

        modelDel2.lb(delIdx_2) = model.lb(delIdx_2);
        modelDel2.ub(delIdx_2) = model.ub(delIdx_2);
        waitbar(iLeth/nLethals, h, [num2str(round(iLeth*100/nLethals)) '% completed...']);
    end
    close(h);
end


function [solutionDel1, solutionDel2, solStatus] = get_serrano_results(modelDel1, modelDel2, p_norm, obj_slack, osenseStr, verbFlag)
% Performs the specified p_norm FBA

    if (nargin < 3 || isempty(p_norm))
        p_norm = 1;
    end

    if (nargin < 4 || isempty(obj_slack))
        obj_slack = 0.05;
    end

    if (nargin < 5 || isempty(osenseStr))
        osenseStr = 'max';
    end

    if (nargin < 6)
        verbFlag = false;
    end

    switch p_norm
        case 0
            solutionDel1 = optimizeCbModel(modelDel1, 'max', 0);
            solutionDel2 = optimizeCbModel(modelDel2, 'max', 0);
        case 1
            solutionDel1 = optimizeCbModel(modelDel1, 'max', 'one');
            solutionDel2 = optimizeCbModel(modelDel2, 'max', 'one');
        case 2
            solutionDel1 = optimizeCbModel(modelDel1, 'max', 2);
            solutionDel2 = optimizeCbModel(modelDel2, 'max', 2);
    end

    if (solutionDel2.stat>0) && (solutionDel1.stat>0)
        solStatus = 1;
    else
        solStatus = -1;
    end
end
