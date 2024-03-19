addpath('core')

% This pipeline prefers using XML files as input.
% So, XML files of the model (with the model ID as the name),
% are expected to be present in the path_to_models location.
model_names = {'e_coli_core', 'iECW_1372', 'iIT341', 'iJO1366', 'iML1515', 'iNJ661', 'iPC815', 'iSSON_1240', 'iYL1228', 'iYS1720', 'STM_v1_0'};
path_to_models = {'examples/e_coli_core/', 'examples/iECW_1372/', 'examples/iIT341/', 'examples/iJO1366/', 'examples/iML1515/', 'examples/iNJ661/', 'examples/iPC815/', 'examples/iSSON_1240/', 'examples/iYL1228/', 'examples/iYS1720/', 'examples/STM_v1_0/'};

%%
fprintf('Reading and converting .xml files to .mat ...\n')
% Read and convert the XML files to MAT files
convert_xml_to_mat(path_to_models, model_names);

%%
% Get the Castle structure with all lethal reactions
% 2 for double lethals and 3 for triple lethals
fprintf('Getting all the Castle objects ...\n')
Castle = getFastSL(model_names, path_to_models, 2);
save('results/Castle.mat', 'Castle');

%%
% Run multiMinRerouting for 0, 1 and 2 norm
fprintf('Performing 0-norm minRerouting ...\n')
Castle0 = multiMinRerouting(Castle, 'zero');
save('results/Castle0.mat', 'Castle0');
fprintf('Performing 1-norm minRerouting ...\n')
Castle1 = multiMinRerouting(Castle, 'one');
save('results/Castle1.mat', 'Castle1');
fprintf('Performing 2-norm minRerouting ...\n')
Castle2 = multiMinRerouting(Castle, 'two');
save('results/Castle2.mat', 'Castle2');

%%
% Final pFBA and PSL-RSL Analysis
pfba_analysis(model_names, path_to_models)
fva_analysis(model_names, path_to_models)
