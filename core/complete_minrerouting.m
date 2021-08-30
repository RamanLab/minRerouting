% addpath('core')

% Recommneded to rename model.mat file to modelNames.mat file
model_names = {'iIT341', 'iEK1008', 'iML1515', 'iNJ661', 'iPC815', 'iYL1228', 'STM_v1_0', 'e_coli_core', 'iJO1366', 'iSSON_1240'};
path_to_models = {'examples/iIT341/', 'examples/iEK1008/', 'examples/iML1515/', 'examples/iNJ661/','examples/iPC815/', 'examples/iYL1228/', 'examples/STM_v1_0/', 'examples/e_coli_core/', 'examples/iJO1366/', 'examples/iSSON_1240/'};

% 2 for double lethals and 3 for triple lethals
Castle = getFastSL(model_names, path_to_models, 2);

% Run multiMinRerouting for 0 norm solution
Castle0 = multiMinRerouting(Castle, 'zero');
% Run multiMinRerouting for 1 norm solution
Castle1 = multiMinRerouting(Castle, 'one');
% Run multiMinRerouting for 2 norm solution
Castle2 = multiMinRerouting(Castle, 'two');

pfba_analysis(model_names, path_to_models)
psl_rsl_analysis(model_names, path_to_models)