function [Castle] = getFastSL(model_names, path_to_models, leth_order)
%getFastSL returns Castle object with multiple models and synthetic lethals 
% Input:
% model_names     List of model IDs e.g. iML1515. Models should be stored
% in path with modelID.mat name
% path_to_models  List of path to directory of model files. Files
% relevant to each model will be stored in this path
% leth_order      Order of highest synthetic lethals 
%
% Output:
% Castle.data   The structure with multiple models and synthetic lethals
%   model       Cobra model structures for each of organism
%   model_name  Names of the models for comparison
%   Jsl         List of single lethals
%   Jdl         List of double lethals
%   Jtl         List of triple lethals
% 
% Omkar Mohite       10 Jul,2018.

disp("Models to be find synthetic lethals for: " + length(model_names))

for j = 1:length(path_to_models)
    model = load(strcat(path_to_models{j}, model_names{j},'.mat'));
    Castle.data(j).model_name = model_names{j};
    Castle.data(j).model = model.(model_names{j});
    clear model;
    
    disp('Eliminating exchange reactions from synthetic lethal search space...')
    eliList= Castle.data(j).model.rxns(find(findExcRxns(Castle.data(j).model)));
    disp(strcat('Identifying synthetic lethals for ', model_names{j}, ' using FastSL... This may take a while'))
    fastSL(Castle.data(j).model,0.01,leth_order,eliList)

    SL_dir = dir('*_Rxn_lethals.mat');
    load(SL_dir.name)
    Castle.data(j).Jsl = Jsl;
    if leth_order>1
        Castle.data(j).Jdl = Jdl;
    end
    if leth_order>2
        Castle.data(j).Jtl = Jtl;
    end   
    disp(strcat('Number synthetic lethals identified for ', model_names{j}, ':'))
    disp(strcat(length(Jsl), ' Single Lethals, ',  length(Jdl), ' Double Lethals')) %Update for triple lethals
    movefile (SL_dir.name,path_to_models{j});
end 
end

