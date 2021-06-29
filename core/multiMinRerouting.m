function [Castle] = multiMinRerouting(Castle, minNorm, obj_slack)
% [Castle] = multiMinRerouting(Castle, minNorm) 
% multiMinRerouting identifies minimal rerouting sets
% for synthetic lethals of multiple models in Castle structure 
% Input:
%   Castle :    The structure with models and synthetic lethals data
%   minNorm :   'two', 'one' or 'zero' identifies minReroutingSets 
%               using minimal L2, L1 or L0 norm formulation of LP
% OPTIONAL
%   obj_slack   Permissible slack on the WT growth rate in the MoMA
%               default value: 0.05

if (nargin < 3 || isempty(obj_slack))
    obj_slack = 0.05;
end

nModels = length(Castle.data);
for j= 1:nModels
    if strcmp(minNorm, 'two')
        disp("Finding set of minimal rerouting reactions using L2 minNorm in " + Castle.data(j).model_name);
        % quadraticMOMA version    
        Castle.data(j).minRe = minReroutingRxns_l2_moma(Castle.data(j).model, Castle.data(j).Jdl, obj_slack);
        fname = "results/" + Castle.data(j).model_name + "/" + "Castle_two_norm.mat";
        data = Castle.data(j);
        save(fname, 'data')
    
    elseif strcmp(minNorm, 'one')
        disp("Finding set of minimal rerouting reactions using L1 minNorm in " + Castle.data(j).model_name);
        % linearMOMA version   
        Castle.data(j).minRe = minReroutingRxns_l1_moma(Castle.data(j).model, Castle.data(j).Jdl, obj_slack);
        fname = "results/" + Castle.data(j).model_name + "/" + "Castle_one_norm.mat";
        data = Castle.data(j);
        save(fname, 'data')
    
    elseif strcmp(minNorm, 'zero')
        disp("Finding set of minimal rerouting reactions using L0 minNorm (iterative formulation) in " + Castle.data(j).model_name);
        % zero-norm LP version   
        Castle.data(j).minRe = minReroutingRxns_l0_moma(Castle.data(j).model, Castle.data(j).Jdl, obj_slack);
        fname = "results/" + Castle.data(j).model_name + "/" + "Castle_zero_norm.mat";
        data = Castle.data(j);
        save(fname, 'data')
    end
end    
end
