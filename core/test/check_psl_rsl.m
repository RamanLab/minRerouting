% model, Jdl
model = load('examples/iJO1366/iJO1366.mat');
model = model.model;
Jdl = readtable('notebooks/cleaned_iJO1366.csv');
solution = optimizeCbModel(model);

a = readtable('results/iJO1366/iJO1366_FVA_one_norm_100.csv')
a(find(~(a.Rxn1_Prod > 0) & ~(a.Rxn2_Prod > 0)), :)

g_max_Fluxes = zeros(numel(Jdl.Rxn_1), 1);
l_max_Fluxes = zeros(numel(Jdl.Rxn_1), 1);
g_min_Fluxes = zeros(numel(Jdl.Rx% model, Jdl
model = load('examples/iJO1366/iJO1366.mat');
model = model.model;
Jdl = readtable('notebooks/cleaned_iJO1366.csv');
solution = optimizeCbModel(model);

a = readtable('results/iJO1366/iJO1366_FVA_one_norm_100.csv')
a(find(~(a.Rxn1_Prod > 0) & ~(a.Rxn2_Prod > 0)), :)

g_max_Fluxes = zeros(numel(Jdl.Rxn_1), 1);
l_max_Fluxes = zeros(numel(Jdl.Rxn_1), 1);
g_min_Fluxes = zeros(numel(Jdl.Rxn_1), 1);
l_min_Fluxes = zeros(numel(Jdl.Rxn_1), 1);

for i = 1:numel(Jdl.Rxn_1)
    rxn1 = Jdl.Rxn_1(i);
    rxn2 = Jdl.Rxn_2(i);

    rxn_id1 = findRxnIDs(model, rxn1);
    rxn_id2 = findRxnIDs(model, rxn2);
    

    % Keep the reaction 1 as not equal to zero and 
    % get the max and min of reaction 2
    
    % Condition 1: rxn1 flux > 0
    temp = struct();
    temp.S = model.S;

    condition = zeros(1, size(model.S, 2));
    condition(rxn_id1) = 1;
    
    osense = [model.csense; 'G'];
    temp.S = [model.S; condition];

    temp.csense = model.csense;
    temp.b = zeros(size(temp.S, 1), 1);

    temp.c = zeros(size(temp.S, 2), 1);
    temp.c(rxn_id2) = 1;

    temp.ub = model.ub;
    temp.lb = model.lb;

    % minimize
    temp.osense = 1;
    minFlux = optimizeCbModel(temp, 'min');
    % maximize
    temp.osense = -1;
    maxFlux = optimizeCbModel(temp, 'max');

    g_min_Fluxes(i) = minFlux.f;
    g_max_Fluxes(i) = maxFlux.f;

    % Condition 2: rxn1 flux < 0
    temp = struct();
    temp.S = model.S;

    condition = zeros(1, size(model.S, 2));
    condition(rxn_id1) = 1;

    osenseStr = [model.csense; 'L'];
    temp.S = [model.S; condition];

    temp.csense = model.csense;
    temp.b = zeros(size(temp.S, 1), 1);

    temp.c = zeros(size(temp.S, 2), 1);
    temp.c(rxn_id2) = 1;

    temp.ub = model.ub;
    temp.lb = model.lb;

    % minimize
    temp.osense = 1;
    minFlux = optimizeCbModel(temp, 'min');
    % maximize
    temp.osense = -1;
    maxFlux = optimizeCbModel(temp, 'max');

    l_min_Fluxes(i) = minFlux.f;
    l_max_Fluxes(i) = maxFlux.f;
end

Jdl.g_min_Fluxes = g_min_Fluxes;
Jdl.g_max_Fluxes = g_max_Fluxes;
Jdl.l_min_Fluxes = l_min_Fluxes;
Jdl.l_max_Fluxes = l_max_Fluxes;

writetable(Jdl, 'psl_test.csv');n_1), 1);
l_min_Fluxes = zeros(numel(Jdl.Rxn_1), 1);

for i = 1:numel(Jdl.Rxn_1)
    rxn1 = Jdl.Rxn_1(i);
    rxn2 = Jdl.Rxn_2(i);

    rxn_id1 = findRxnIDs(model, rxn1);
    rxn_id2 = findRxnIDs(model, rxn2);
    

    % Keep the reaction 1 as not equal to zero and 
    % get the max and min of reaction 2
    
    % Condition 1: rxn1 flux > 0
    temp = struct();
    temp.S = model.S;

    condition = zeros(1, size(model.S, 2));
    condition(rxn_id1) = 1;
    
    osense = [model.csense; 'G'];
    temp.S = [model.S; condition];

    temp.csense = model.csense;
    temp.b = zeros(size(temp.S, 1), 1);

    temp.c = zeros(size(temp.S, 2), 1);
    temp.c(rxn_id2) = 1;

    temp.ub = model.ub;
    temp.lb = model.lb;

    % minimize
    temp.osense = 1;
    minFlux = optimizeCbModel(temp, 'min');
    % maximize
    temp.osense = -1;
    maxFlux = optimizeCbModel(temp, 'max');

    g_min_Fluxes(i) = minFlux.f;
    g_max_Fluxes(i) = maxFlux.f;

    % Condition 2: rxn1 flux < 0
    temp = struct();
    temp.S = model.S;

    condition = zeros(1, size(model.S, 2));
    condition(rxn_id1) = 1;

    osenseStr = [model.csense; 'L'];
    temp.S = [model.S; condition];

    temp.csense = model.csense;
    temp.b = zeros(size(temp.S, 1), 1);

    temp.c = zeros(size(temp.S, 2), 1);
    temp.c(rxn_id2) = 1;

    temp.ub = model.ub;
    temp.lb = model.lb;

    % minimize
    temp.osense = 1;
    minFlux = optimizeCbModel(temp, 'min');
    % maximize
    temp.osense = -1;
    maxFlux = optimizeCbModel(temp, 'max');

    l_min_Fluxes(i) = minFlux.f;
    l_max_Fluxes(i) = maxFlux.f;
end

Jdl.g_min_Fluxes = g_min_Fluxes;
Jdl.g_max_Fluxes = g_max_Fluxes;
Jdl.l_min_Fluxes = l_min_Fluxes;
Jdl.l_max_Fluxes = l_max_Fluxes;

writetable(Jdl, 'psl_test.csv');