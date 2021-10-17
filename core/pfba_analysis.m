function pfba_analysis(model_names, path_to_models)
% pfba_analysis performs a pFBA on the models provided.
% The double lethals reactions and their classes are saved 
% in the form of .csv files in their respective folders.
% 
% Input:
% model_names     List of model IDs e.g. iML1515. Models should be stored
%                 in path with modelID.mat name
% path_to_models  List of path to directory of model files. Files
%                 relevant to each model will be stored in this path
%

for j = 1:numel(model_names)
    fprintf('Performing pFBA analysis for %s ...\n', model_names{j});
    model = load(strcat(path_to_models{j}, model_names{j},'.mat'));
    model = model.model;
    Castle.data = load(strcat(path_to_models{j}, model_names{j},'_Rxn_lethals.mat'));

    % Performing pFBA with geneoption as 0
    % 0 ensures that flux is reduced through all reactions
    % Default: 1 (only gene-associated fluxes)
    [~, RxnClasses, ~] = pFBA(model, 'geneoption', 0);
    try
        save(strcat('results/', model_names{j}, '/pfba/RxnClasses.mat'), 'RxnClasses')
    catch
        try
            mkdir(strcat('results/', model_names{j}, '/pfba/'))
        catch 
        end
        save(strcat('results/', model_names{j}, '/pfba/RxnClasses.mat'), 'RxnClasses')
    end

    % Save the list of all unique reactions in the Jdl
    % containers.Map is used to map the reactions to their
    % respective reaction classes (similar to dict in Python)
    pfba_rxn_analysis = Castle.data.Jdl;
    unique_rxns = unique(pfba_rxn_analysis);
    unique_rxn_mapping = containers.Map();

    % Map the reactions to their classes
    for i = 1:numel(unique_rxns)
        if sum(ismember(RxnClasses.Essential_Rxns, unique_rxns(i)))
            unique_rxn_mapping(unique_rxns{i}) = 'Essential_Rxns';
        elseif sum(ismember(RxnClasses.pFBAOpt_Rxns, unique_rxns(i)))
            unique_rxn_mapping(unique_rxns{i}) = 'pFBAOpt_Rxns';
        elseif sum(ismember(RxnClasses.ELE_Rxns, unique_rxns(i)))
            unique_rxn_mapping(unique_rxns{i}) = 'ELE_Rxns';
        elseif sum(ismember(RxnClasses.MLE_Rxns, unique_rxns(i)))
            unique_rxn_mapping(unique_rxns{i}) = 'MLE_Rxns';
        elseif sum(ismember(RxnClasses.ZeroFlux_Rxns, unique_rxns(i)))
            unique_rxn_mapping(unique_rxns{i}) = 'ZeroFlux_Rxns';
        elseif sum(ismember(RxnClasses.Blocked_Rxns, unique_rxns(i)))
            unique_rxn_mapping(unique_rxns{i}) = 'Blocked_Rxns';
        end
    end

    % Adding the classes to the cell matrix
    for i = 1:size(pfba_rxn_analysis)
        pfba_rxn_analysis{i,3} = unique_rxn_mapping(pfba_rxn_analysis{i,1});
        pfba_rxn_analysis{i,4} = unique_rxn_mapping(pfba_rxn_analysis{i,2});
    end
    clear i

    % Converting the cell matrix to table and renaming columns
    pfba_table = cell2table(pfba_rxn_analysis, 'VariableNames',{'Rxn_1' 'Rxn_2' 'Rxn_1_Class' 'Rxn_2_Class'});
    
    % Saving the table in .csv format
    writetable(pfba_table, strcat(model_names{j}, '_pFBA.csv'));
    % Move the file back to the original path of the model
    save_path = strcat('results/', model_names{j});
    movefile(strcat(model_names{j}, '_pFBA.csv'), save_path);
end
end
