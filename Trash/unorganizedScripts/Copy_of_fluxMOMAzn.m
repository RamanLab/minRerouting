function [solutionWT] = fluxMOMAzn( modeldel_1,modeldel_2 )
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

[nMets,nRxns] = size(modeldel_1.S);


% Define the constraints structure
    constraint.A = [modeldel_1.S;modeldel_2.S];
    constraint.b = [modeldel_1.b;modeldel_2.b];
    constraint.csense(1:2*nMets,1) = 'E'; 
    constraint.lb = [modeldel_1.lb;modeldel_2.lb];
    constraint.ub = [modeldel_1.ub;modeldel_2.ub];
    
    params.epsilon=eps;
    
    % Call the sparse LP solver
    solutionL0 = sparseLP('all',constraint);
    
    %Store results
    solution.stat   = solutionL0.stat;
    solution.full   = solutionL0.x;
    solution.dual   = [];
    solution.rcost  = [];        
%     solution.dual=solution.dual(1:m,1);
    solutionWT.x = solution.full(1:nRxns)+fluxDel;
end

