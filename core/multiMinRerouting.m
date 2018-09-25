function [Castle] = multiMinRerouting(Castle, minNorm)
% [Castle] = multiMinRerouting(Castle, minNorm) 
% multiMinRerouting identifies minimal rerouting sets for synthetic lethals of multiple models in Castle structure 
%   Input:
%   Castle : The structure with models and synthetic lethals data
%   minNorm : 'two', 'one' or 'zero' identifies minReroutingSets using minimal L2, L1 or L0 norm formulation of LP

nModels = length(Castle.data);
for j= 1:nModels
    if strcmp(minNorm, 'two')
        disp("Finding set of minimal rerouting reactions using L2 minNorm in " + Castle.data(j).model_name);
        Castle.data(j).minRe = minReroutingRxns_l2_moma(Castle.data(j).model,Castle.data(j).Jdl);    %quadraticMOMA version    
    elseif strcmp(minNorm, 'one')
        disp("Finding set of minimal rerouting reactions using L1 minNorm in " + Castle.data(j).model_name);
        Castle.data(j).minRe = minReroutingRxns_l1_moma(Castle.data(j).model,Castle.data(j).Jdl);    %linearMOMA version   
    elseif strcmp(minNorm, 'zero')
        disp("Finding set of minimal rerouting reactions using L0 minNorm (iterative formulation) in " + Castle.data(j).model_name);
        Castle.data(j).minRe = minReroutingRxns_l0_moma(Castle.data(j).model,Castle.data(j).Jdl); % Iterative formulation    
    end
end    
end

