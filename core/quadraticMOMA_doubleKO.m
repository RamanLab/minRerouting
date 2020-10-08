function [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = quadraticMOMA_doubleKO(modelDel1, modelDel2, osenseStr, minFluxFlag, verbFlag)
% Performs a quadratic version of the MOMA (minimization of metabolic
% adjustment) approach upgraded for comparing double knockouts
%
% USAGE:
%
%    [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = quadraticMOMA_doubleKO(modelDel1, modelDel2, osenseStr, minFluxFlag, verbFlab)
%
% INPUTS:
%    modelDe11:         Deletion strain model for Rxn1
%    modelDel2:         Deletion strain model for Rxn2
%
% OPTIONAL INPUTS:
%    osenseStr:        Maximize ('max') / minimize ('min') (Default = 'max')
%    minFluxFlag:      Minimize the absolute value of fluxes in the optimal MOMA
%                      solution (Default = false)
%    verbFlag:         Verbose output (Default = false)
%
% OUTPUTS:
%    solutionDel1:     Deletion 1 solution structure
%    solutionDel2:     Deletion 2 solution structure
%    totalFluxDiff:    Value of the linear MOMA objective, i.e. :math:`\sum |v_{del1}-v_{del2}|`
%    solStatus:        Solution status - solves the problem: (`f_wt` is the optimal wild type objective value found by FBA)
%
% .. math::
%     min ~&~  \sum |v_{del1} - v_{del2}|_2 \\
%         ~&~ S_{del1}v_{del1} = 0 \\
%         ~&~ lb_{del1} \leq v_{del1} \leq ub_{del1} \\
%         ~&~ c_{del1}^T v_{del1} = f_{del1} \\
%         ~&~ S_{del2}v_{del2} = 0 \\
%         ~&~ lb_{del2} \leq v_{del2} \leq ub_{del2}
%         ~&~ c_{del2}^T v_{del2} = f_{del2} \\
%
% NOTE:
%
%    1) This formulation allows for selecting the most appropriate
%    optimal FBA solutions for both knockout strains as the starting point as opposed to
%    picking an arbitrary starting point (original MOMA implementation).
%
%    2) The reaction sets in the two models do not have to be equal as long as
%    there is at least one reaction in common
%
% .. Author: - Markus Herrgard 11/7/06
% ..           Omkar Satyavan Mohite 06-08-2018

if (nargin < 3 || isempty(osenseStr))
    osenseStr = 'max';
end
if (nargin < 4 || isempty(minFluxFlag))
    minFluxFlag = false;
end
if (nargin < 5)
    verbFlag = false;
end

% LP solution tolerance
global CBT_LP_PARAMS
if (exist('CBT_LP_PARAMS', 'var'))
    if isfield(CBT_LP_PARAMS, 'objTol')
        tol = CBT_LP_PARAMS.objTol;
    else
        tol = 1e-6;
    end
else
    tol = 1e-6;
end

[nMets1,nRxns1] = size(modelDel1.S);
[nMets2,nRxns2] = size(modelDel2.S);

% Match model reaction sets
commonRxns = ismember(modelDel1.rxns,modelDel2.rxns);
nCommon = sum(commonRxns);
if (nCommon == 0)
    error('No common rxns in the models');
end

solutionDel1.f = [];
solutionDel1.x = [];
solutionDel1.stat = -1;
solutionDel2.f = [];
solutionDel2.x = [];
solutionDel2.stat = -1;

if (verbFlag)
    fprintf('Solving FBA for deletion 1 strain: %d constraints %d variables ', nMets1, nRxns1);
end
% Solve wt problem
solutionDel1 = optimizeCbModel(modelDel1,osenseStr, 'one');

if (verbFlag)
    fprintf('%f seconds\n',solutionDel1.time);
end
% Round off solution to avoid numerical problems
if (strcmp(osenseStr,'max'))
    objValDel1 = floor(solutionDel1.f/tol)*tol;
else
    objValDel1 = ceil(solutionDel1.f/tol)*tol;
end

if (verbFlag)
    fprintf('Solving FBA for deletion 2 strain: %d constraints %d variables ',nMets1,nRxns2);
end
% Solve wt problem
solutionDel2 = optimizeCbModel(modelDel2,osenseStr, 'one');

if (verbFlag)
    fprintf('%f seconds\n',solutionDel2.time);
end
% Round off solution to avoid numerical problems
if (strcmp(osenseStr,'max'))
    objValDel2 = floor(solutionDel2.f/tol)*tol;
else
    objValDel2 = ceil(solutionDel2.f/tol)*tol;
end

% Variables in the following problem are
% x = [v1;v2;delta]
% where v1 = wild type flux vector
%       v2 = deletion strain flux vector
%       delta = v1 - v2

if (solutionDel1.stat > 0 && solutionDel2.stat > 0)
    % Construct the LHS matrix
    % Rows:
    % 1: Sdel1*v1 = 0 for the deletion strain 1
    % 2: Sdel2*v2 = 0 for the deletion strain 2
    % 3: delta+ >= v1-v2
    % 4: delta- >= v2-v1
    % 5: c'v1 >= 0.9*f1 (deletion strain 1) (10 % slack on obj)
    % 6: c'v2 >= 0.9*f2 (deletion strain 2)
    obj_slack = 0.1;
    deltaMat = createDeltaMatchMatrix(modelDel1.rxns,modelDel2.rxns);
    deltaMat = deltaMat(1:nCommon,1:(nRxns1+nRxns2+nCommon));
    A = [modelDel1.S sparse(nMets1,nRxns2+nCommon);
        sparse(nMets2,nRxns1) modelDel2.S sparse(nMets2,nCommon);
        deltaMat;
        modelDel1.c' sparse(1,nRxns2+nCommon);
        sparse(1,nRxns1) modelDel2.c' sparse(1,nCommon)];
    
    % Construct the RHS vector
    b = [zeros(nMets1+nMets2+nCommon,1);(1-obj_slack)*objValDel1 ;(1-obj_slack)*objValDel2];
    
    c = [zeros(nRxns1+nRxns2+nCommon,1)];
    
    % Construct the ub/lb
    % delta [-10000 10000]
    lb = [modelDel1.lb;modelDel2.lb;-10000*ones(nCommon,1)];
    ub = [modelDel1.ub;modelDel2.ub;10000*ones(nCommon,1)];
    
    % Construct the constraint direction vector (G for delta's, E for
    % everything else)
    csense(1:(nMets1+nMets2+nCommon)) = 'E';
    if (strcmp(osenseStr,'max'))
        csense(end+1) = 'G';
        csense(end+1) = 'G';
    else
        csense(end+1) = 'L';
        csense(end+1) = 'L';
    end
    
    % F matrix
    F = [sparse(nRxns1+nRxns2,nRxns1+nRxns2+nCommon);
        sparse(nCommon,nRxns1+nRxns2) 2*eye(nCommon)];
    
    if (verbFlag)
        fprintf('Solving quadratic MOMA upgraded for double knockouts: %d constraints %d variables ',size(A,1),size(A,2));
    end
    
    % Solve the quadraticMOMA problem
    [QPproblem.A,QPproblem.b,QPproblem.F,QPproblem.c,QPproblem.lb,QPproblem.ub,QPproblem.csense,QPproblem.osense] = deal(A,b,F,c,lb,ub,csense,1);
    %QPsolution = solveCobraQP(QPproblem,[],verbFlag-1);
    QPsolution = solveCobraQP(QPproblem);
    
    if (verbFlag)
        fprintf('%f seconds\n',QPsolution.time);
    end
    
    % Get the solution(s)
    if QPsolution.stat == 1
        solutionDel1.x = QPsolution.full(1:nRxns1);
        solutionDel2.x = QPsolution.full((nRxns1+1):(nRxns1+nRxns2));
        
        solutionDel1.f = sum(modelDel2.c.*solutionDel1.x);
        solutionDel2.f = sum(modelDel2.c.*solutionDel2.x);
        
        totalFluxDiff = sum((solutionDel1.x-solutionDel2.x).^2);
    else
        totalFluxDiff = [];
    end
    solutionDel1.stat = QPsolution.stat;
    solutionDel2.stat = QPsolution.stat;
    solStatus = QPsolution.stat;
    
elseif solutionDel1.stat <= 0
    warning('Deletion 1 strain FBA problem is infeasible or unconstrained');
    solStatus = solutionDel1.stat;
    totalFluxDiff = [];
else
    warning('Deletion 2 strain FBA problem is infeasible or unconstrained');
    solStatus = solutionDel2.stat;
    totalFluxDiff = [];
end
end
