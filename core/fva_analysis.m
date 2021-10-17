function fva_analysis(model_names, path_to_models)
% fva_analysis performs a FVA on the models provided.
% The flux variability of the double lethal pairs and a zero 
% norm solution of the model is saved in the form of
% .csv files in their respective folders.
% 
% Input:
% model_names     List of model IDs e.g. iML1515. Models should be stored
%                 in path with modelID.mat name
% path_to_models  List of path to directory of model files. Files
%                 relevant to each model will be stored in this path
%

for i = 1:numel(model_names)
    % Get the lethal reactions data
    data = load(strcat(path_to_models{i}, model_names{i}, '_Rxn_lethals.mat'));
    model = load(strcat(path_to_models{i}, model_names{i}, '.mat'));
    model = model.model;

    % Perform an FVA on the model
    fprintf('Finding Flux Variability of %s ... ', model_names{i})
    
    % Exclude exchange reactions
    rxnList = model.rxns(~findExcRxns(model));
    % Increase the level of verbosity and make the optimization as 1-norm
    [minFlux, maxFlux] = fluxVariability(model, 100, 'max', rxnList, 1, 1, '1-norm');
    
    % Get indices of the first reactions
    [~,positions1] = ismember(data.Jdl(:,1), rxnList);
    maxFlux1 = maxFlux(positions1);
    minFlux1 = minFlux(positions1);

    % Get indices of the second reactions
    [~,positions2] = ismember(data.Jdl(:,2), rxnList);
    maxFlux2 = maxFlux(positions2);
    minFlux2 = minFlux(positions2);

    % Convert the flux values to a tabular form
    fvaResultAnalysis = data.Jdl;
    fvaResultAnalysis = [fvaResultAnalysis, array2table(minFlux1)];
    fvaResultAnalysis = [fvaResultAnalysis, array2table(maxFlux1)];
    fvaResultAnalysis = [fvaResultAnalysis, array2table(minFlux2)];
    fvaResultAnalysis = [fvaResultAnalysis, array2table(maxFlux2)];
    fvaResultAnalysis.Properties.VariableNames = {'Rxn_1' 'Rxn_2' 'Rxn_1_Min' 'Rxn_1_Max' 'Rxn_2_Min' 'Rxn_2_Max'};
    fvaResultAnalysis = fvaResultAnalysis(:,[1,3,4,2,5,6]);
    % Get an one norm FBA solution for the model
    sol = optimizeCbModel(model, 'max', 'one');
    v1 = sol.v(positions1);
    v2 = sol.v(positions2);
    
    fvaResultAnalysis = [fvaResultAnalysis, array2table(v1)];
    fvaResultAnalysis = [fvaResultAnalysis, array2table(v2)];
    
    % Get the product too
    fvaResultAnalysis.Rxn1_Prod = maxFlux1.*minFlux1;
    fvaResultAnalysis.Rxn2_Prod = maxFlux2.*minFlux2;
    
    % Save the table in .csv format
    writetable(fvaResultAnalysis, strcat(model_names{i}, '_FVA_one_norm_100.csv'));
    
    % Move the file back to the original path of the model
    save_path = strcat('results/', model_names{i});
    movefile(strcat(model_names{i}, '_FVA_one_norm_100.csv'), save_path);
    fprintf('Done!\n')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Analysis Part 2; Case where the category of 
    % double lethals is ambiguous.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    solution = optimizeCbModel(model);
    v_bio_opt = solution.f;
    Jdl = fvaResultAnalysis(find(~((fvaResultAnalysis.Rxn1_Prod > 0) & (fvaResultAnalysis.Rxn2_Prod > 0))), :);

    g_max_Fluxes = zeros(numel(Jdl.Rxn_1), 1);
    l_max_Fluxes = zeros(numel(Jdl.Rxn_1), 1);
    g_min_Fluxes = zeros(numel(Jdl.Rxn_1), 1);
    l_min_Fluxes = zeros(numel(Jdl.Rxn_1), 1);

    for j = 1:numel(Jdl.Rxn_1)
        rxn1 = Jdl.Rxn_1(j);
        rxn2 = Jdl.Rxn_2(j);

        rxn_id1 = findRxnIDs(model, rxn1);
        rxn_id2 = findRxnIDs(model, rxn2);
        
        % Keep the reaction 1 as not equal to zero and 
        % get the max and min of reaction 2
        % The proble solved by optimizesCbModel is
        %  max/min  c^T v
        % s.t.  Sv = b
        %       lb <= v <= ub
        
        % Condition 1: rxn1 flux > 0
        temp = struct();
        temp.S = model.S;

        condition = zeros(1, size(model.S, 2));
        condition(rxn_id1) = 1;
        temp.S = [model.S; condition];
        
        csense = [model.csense; 'G'];
        temp.csense = csense;
        temp.b = zeros(size(temp.S, 1), 1);

        temp.c = zeros(size(temp.S, 2), 1);
        temp.c(rxn_id2) = 1;

        temp.ub = model.ub;
        temp.lb = model.lb;

        % Ensure that the vbio is same as WT        
        temp.lb(find(model.c)) = solution.f;
        temp.ub(find(model.c)) = solution.f;

        % minimize
        temp.osense = 1;
        minFlux = optimizeCbModel(temp, 'min');
        % maximize
        temp.osense = -1;
        maxFlux = optimizeCbModel(temp, 'max');

        g_min_Fluxes(j) = minFlux.f;
        g_max_Fluxes(j) = maxFlux.f;

        % Condition 2: rxn1 flux < 0
        temp = struct();
        temp.S = model.S;

        condition = zeros(1, size(model.S, 2));
        condition(rxn_id1) = 1;
        temp.S = [model.S; condition];

        csense = [model.csense; 'L'];
        temp.csense = csense;
        temp.b = zeros(size(temp.S, 1), 1);

        temp.c = zeros(size(temp.S, 2), 1);
        temp.c(rxn_id2) = 1;

        temp.ub = model.ub;
        temp.lb = model.lb;
        
        % Ensure that the vbio is same as WT        
        temp.lb(find(model.c)) = solution.f;
        temp.ub(find(model.c)) = solution.f;

        % minimize
        temp.osense = 1;
        minFlux = optimizeCbModel(temp, 'min');
        % maximize
        temp.osense = -1;
        maxFlux = optimizeCbModel(temp, 'max');

        l_min_Fluxes(j) = minFlux.f;
        l_max_Fluxes(j) = maxFlux.f;
    end

    Jdl.g_min_Fluxes = g_min_Fluxes;
    Jdl.g_max_Fluxes = g_max_Fluxes;
    Jdl.l_min_Fluxes = l_min_Fluxes;
    Jdl.l_max_Fluxes = l_max_Fluxes;

    writetable(Jdl, strcat('results/', model_names{i}, '/', model_names{i}, '_ambiguous_FVA_100.csv'));
end
end
