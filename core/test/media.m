model_names = {'iJO1366'};
path_to_models = {'examples/iJO1366/'};

for j = 1:numel(model_names)
    model = load(strcat(path_to_models{j}, model_names{j},'.mat'));
    model = model.model;
    model.rxns(model.lb~=0  & model.lb~=-1000)
    index = findRxnIDs(model, 'EX_glc__D_e');
    atpm_rxn_id = findRxnIDs(model, 'ATPM');
    atpm_lb = model.lb(atpm_rxn_id);
    atpm_ub = model.ub(atpm_rxn_id);
    
    model.lb(index) = -10;
    model.ub(index) = 0;
    save(strcat(path_to_models{j}, model_names{j},'_minimal.mat'), 'model');
end

model_names = {'iJO1366_minimal'};
path_to_models = {'examples/iJO1366/'};
Castle = getFastSL(model_names, path_to_models, 2);

Castle0 = multiMinRerouting(Castle, 'zero');
% save('results/Castle0.mat', 'Castle0')
% Run multiMinRerouting for 1 norm solution
Castle1 = multiMinRerouting(Castle, 'one');
% save('results/Castle1.mat', 'Castle1')
% Run multiMinRerouting for 2 norm solution
Castle2 = multiMinRerouting(Castle, 'two');
% save('results/Castle2.mat', 'Castle2')


model_names = {'iJO1366'};
path_to_models = {'examples/iJO1366/'};
pfba_analysis(model_names, path_to_models)
fva_analysis(model_names, path_to_models)
