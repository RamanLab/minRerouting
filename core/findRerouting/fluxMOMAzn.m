function [solutionWT] = fluxMOMAzn( model,fluxDel )
% fluxMOMA performs zero norm version of MOMA for a given model and flux
% vector
%% [solutionWT] = MOMA(model,fluxDel)
%
% INPUT
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         Reaction Names
% fluxDel          Deletion strain flux vector
%
%OUTPUTS
% solution        Wild-type solution structure

% MOMA to find closest flux to given vector using zero norm LP
%               min ||V'||_0
%               s.t. S*V' = +/-SV = 0  
%               lb-fluxDel <= V' <= ub-fluxDel
%               V'=V-fluXDel   
%
% Omkar 15/03/2017

zeroNormApprox = 'cappedL1'; %This is default apporximation used in optimizeCbModel

[nMets,nRxns] = size(model.S);


% Define the constraints structure
    constraint.A = [model.S];
    constraint.b = [model.b];
    constraint.csense(1:nMets,1) = 'E'; 
    constraint.lb = model.lb-fluxDel;
    constraint.ub = model.ub-fluxDel;
    
%     params.epsilon=10e-4;
    
    % Call the sparse LP solver
    solutionL0 = sparseLP('cappedL1',constraint);
    
    %Store results
    solution.stat   = solutionL0.stat;
    solution.full   = solutionL0.x;
    solution.dual   = [];
    solution.rcost  = [];        
%     solution.dual=solution.dual(1:m,1);
    solutionWT.x = solution.full(1:nRxns)+fluxDel;
end

