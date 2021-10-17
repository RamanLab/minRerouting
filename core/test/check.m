model_names = {'iJO1366'};
path_to_models = {'examples/iJO1366/'};

for j = 1:numel(model_names)
    model = load(strcat(path_to_models{j}, model_names{j},'.mat'));
    model = model.model;
    glu_rxn_id = findRxnIDs(model, 'EX_glc__D_e');
    atpm_rxn_id = findRxnIDs(model, 'ATPM');
    
    try
        atpm_lb = model.lb(atpm_rxn_id);
        atpm_ub = model.ub(atpm_rxn_id);
    catch
        % do nothing
    end

    indices = find(model.lb~=0  & model.lb~=-1000);
    model.lb(indices) = 0;
    
    model.lb(glu_rxn_id) = -10;
    model.ub(glu_rxn_id) = 0;
    model.lb(atpm_rxn_id) = atpm_lb;
    model.ub(atpm_rxn_id) = atpm_ub;

   pfba_analysis(model_names, path_to_models)
   
end