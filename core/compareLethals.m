function [Castle] = compareLethals(Castle, leth_order)
% compareDL compares synthetic lethals across multiple models  
% Lethality status : Status of presence and lethality is assigned to each 
% reactions in respective model
% S: Single Lethal,D: Double Lethal, T: Triple Lethal, N: Not part of 
% synthetic lethal (S,D or T), A: Absent
% For example: If reaction 'rxnA' is single lethal in model1 and part of 
% double in model2: String of 'S, D' is assigned to rxnA row. 
%
% Input : 
% Castle.data   The structure with multiple models and synthetic lethals
%   model       Cobra model structures for each of organism
%   model_name  Names of the models for comparison
%   Jsl         List of single lethals
%   Jdl         List of double lethals
%   Jtl         List of triple lethals
% leth_order   Order of highest synthetic lethals 
%
% Output :
% Castle      The structure with multiple models and synthetic lethals   
%   data        Field same as input   
%   lethStatus  mat file with status of lethality for union of reactions in
%   each model
% 
%
% Omkar Mohite       10 Jul,2018.

nModels = length(Castle.data);
disp("Models to be compared for synthetic lethals: " + nModels)

% Find unique reactions in all models
Uniq_Rxns = [];
for j = 1:nModels
    Uniq_Rxns = [Uniq_Rxns; Castle.data(j).model.rxns];
end 
Uniq_Rxns = unique(Uniq_Rxns);
disp("Union of reactions in all the models: " + length(Uniq_Rxns))
Castle.all_rxns = Uniq_Rxns;
disp("Finding the status of lethality for every reaction in each of the models ...")

lethStatus = cell(length(Uniq_Rxns), nModels);
lethCount = zeros(length(Uniq_Rxns), 4);

for i = 1:length(Uniq_Rxns)
    S_cnt = 0;
    D_cnt = 0;
    N_cnt = 0;
    A_cnt = 0;
    for j = 1:nModels
        if ismember(Uniq_Rxns(i), Castle.data(j).Jsl)
            lethStatus{i,j} = 'S';
            S_cnt = S_cnt + 1;
        elseif ismember(Uniq_Rxns(i), Castle.data(j).Jdl)
            lethStatus{i,j} = 'D';
            D_cnt = D_cnt + 1;
        elseif leth_order == 2 && ismember(Uniq_Rxns(i), Castle.data(j).model.rxns)
            lethStatus{i,j} = 'N';
            N_cnt = N_cnt + 1;
        % elseif leth_order == 3 && ismember(Uniq_Rxns(i),Castle.data(j).Jtl)
        %     lethStatus{i,j}='T';
        % elseif leth_order == 3 && ismember(Uniq_Rxns(i),Castle.data(j).model.rxns)
        %     lethStatus{i,j}='N';
        else
            lethStatus{i,j} = 'A';
            A_cnt = A_cnt + 1;
        end
    end
    lethCount(i,:) = [S_cnt, D_cnt, N_cnt, A_cnt];
end
Castle.lethStatus = lethStatus;
Castle.lethCount = lethCount;

Castle.core.Single = Uniq_Rxns(find(lethCount(:,1) == 6));
Castle.core.Double = Uniq_Rxns(find(lethCount(:,2) == 6));
Castle.core.Nonessential = Uniq_Rxns(find(lethCount(:,3) == 6));
tmp = Uniq_Rxns(find(lethCount(:,1)+lethCount(:,2) == 6));
tmp = tmp(~ismember(tmp, Castle.core.Single));
Castle.core.Single_to_Double = tmp(~ismember(tmp, Castle.core.Double));

% Similarity in synthetic lethals across models
Castle.simTable.SingleCnt = zeros(nModels, nModels);
for i = 1:nModels
    for j = 1:nModels
        Castle.simTable.SingleCnt(i,j) = length(intersect(Castle.data(i).Jsl, Castle.data(j).Jsl));
        Castle.simTable.Single(i,j).rxns = intersect(Castle.data(i).Jsl, Castle.data(j).Jsl);
    end
end

% Sinlge lethal similarity table of rxn



% Pairwise comparison of models
% New function

end
