function convert_xml_to_mat(path_to_models, model_names)
%% convert_xml_to_mat(path_to_models, model_names)
% INPUT
%   - path_to_models: Path from the parent folder to 
%                     where the model XML files are stored
%   - model_names: Names of the models (XML file names)
% RETURNS
%   - The corresponding model mat files are saved in the 
%     respective folders (path folders)
% 
% N Sowmya Manojna      25/09/2021

    for j = 1:length(path_to_models)
        model = readCbModel(strcat(path_to_models{j}, model_names{j},'.xml'));
        fprintf('Model: %s\n', model_names{j})
        strcat(path_to_models{j}, model_names{j}, '.mat');
        save(strcat(path_to_models{j}, model_names{j}, '.mat'), 'model');
    end
end