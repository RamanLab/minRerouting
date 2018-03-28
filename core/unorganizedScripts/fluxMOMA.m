function [solutionWT, totalFluxDiff] = fluxMOMA(modelWT,fluxDel)
%MOMA Performs a quadratic version of the MOMA (minimization of
%metabolic adjustment) approach 
% Improved MOMA for flux input of del strain 
%
% [solutionDel,solutionWT,totalFluxDiff] = MOMA(modelWT,fluxDel)
%
%INPUTS
% modelWT          Wild type model
% fluxDel          Deletion strain flux vector
%
%OUTPUTS
% solutionWT        Wild-type solution structure
% totalFluxDiff     Value of the linear MOMA objective, i.e.
%                   sum(flux_del-v_wt)^2

% MOMA that avoids problems with alternative optima 
%
%    min sum(flux_del - v_wt)^2
%     S_wt*v_wt = 0
%     lb_wt <= v_wt <= ub_wt
%
% Omkar 03/01/2017

[nMets,nRxns] = size(modelWT.S);

solutionWT.f = [];
solutionWT.x = [];

% Variables in the following problem are
% x = [v1;delta]
% where v1 = wild type flux vector
%       v2 = deletion strain flux vector defined
%       delta = v1 - v2

        b = zeros(nMets,1);
        A = modelWT.S;
        c = -1*fluxDel;
        F = eye(nRxns);
        lb = modelWT.lb;
        ub = modelWT.ub;
        csense(1:nMets) = 'E';

     
    % Solve the linearMOMA problem
    [QPproblem.A,QPproblem.b,QPproblem.F,QPproblem.c,QPproblem.lb,QPproblem.ub,QPproblem.csense,QPproblem.osense] = deal(A,b,F,c,lb,ub,csense,1);
    %QPsolution = solveCobraQP(QPproblem,[],verbFlag-1);
    QPsolution = solveCobraQP(QPproblem);

    % Get the solution(s)
        solutionWT.x = QPsolution.full;
        solutionWT.f = sum(modelWT.c.*solutionWT.x);
        totalFluxDiff = sum((solutionWT.x-fluxDel).^2);
        
end    
    