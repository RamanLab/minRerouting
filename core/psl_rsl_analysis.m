function psl_rsl_analysis(model_names, path_to_models)
% psl_rsl_analysis performs a FVA on the models provided.
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
    model = model.(model_names{i});

    % Perform an FVA on the model
    fprintf('Finding Flux Variability of %s ...\n', model_names{i})
    [minFlux, maxFlux] = fluxVariability(model, [], [], [], 0, 1, '1-norm');

    % Get indices of the first reactions
    positions1 = findRxnIDs(model, data.Jdl(:,1));
    maxFlux1 = maxFlux(positions1);
    minFlux1 = minFlux(positions1);

    % Get indices of the second reactions
    positions2 = findRxnIDs(model, data.Jdl(:,2));
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
    
    % Save the table in .csv format
    writetable(fvaResultAnalysis, strcat(model_names{i}, '_FVA_one_norm.csv'));
    % Move the file back to the original path of the model
    movefile(strcat(model_names{i}, '_FVA_one_norm.csv'), path_to_models{i});
    fprintf('Done!\n')
end
end
