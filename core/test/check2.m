model_names = {'e_coli_core', 'iIT341', 'iJO1366', 'iML1515', 'iNJ661', 'iPC815', 'iYL1228', 'STM_v1_0'};
path_to_models = {'examples/e_coli_core/', 'examples/iIT341/', 'examples/iJO1366/', 'examples/iML1515/', 'examples/iNJ661/', 'examples/iPC815/', 'examples/iYL1228/', 'examples/STM_v1_0/'};

for j = 1:numel(model_names)
    model = load(strcat(path_to_models{j}, model_names{j},'.mat'));
    model = model.model;

    model.modelID
    a = optimizeCbModel(model);
    a.f

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

    try
        model.lb(atpm_rxn_id) = atpm_lb;
        model.ub(atpm_rxn_id) = atpm_ub;
    catch
    end

    b = optimizeCbModel(model);
    b.f
end