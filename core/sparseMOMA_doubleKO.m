function [solutionDel1, solutionDel2, solStatus] = sparseMOMA_doubleKO(modelDel1, modelDel2, osenseStr, minFluxFlag, verbFlag)
% Performs a sparse version of the MOMA (minimization of metabolic
% adjustment) approach upgraded for comparing double knockouts
%
% USAGE:
%
%    [solutionDel1, solutionDel2, totalFluxDiff, solStatus] = sparseMOMA_doubleKO(modelDel1, modelDel2, osenseStr, minFluxFlag, verbFlab)
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
%     min ~&~  \sum |v_{del1} - v_{del2}|_0 \\
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

if (nargin <3 || isempty(osenseStr))
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
    fprintf('Solving FBA for deletion 1 strain: %d constraints %d variables ',nMets1,nRxns1);
end
% Solve wt problem
solutionDel1 = optimizeCbModel(modelDel1,osenseStr, 'zero');

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
solutionDel2 = optimizeCbModel(modelDel2,osenseStr, 'zero');

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
% x = [v1;v2;delta+;delta-]
% where v1 = deletion strain 1 flux vector
%       v2 = deletion strain 2 flux vector
%       delta+ = v1 - v2
%       delta- = v2 - v1
% First solve LP using linearLP and use these as constraint to find minimal set of reactions

if (solutionDel1.stat > 0 && solutionDel2.stat > 0)
    % Construct the LHS matrix
    % Rows:
    % 1: Sdel1*v1 = 0 for the deletion strain 1
    % 2: Sdel2*v2 = 0 for the deletion strain 2
    % 3: delta+ >= v1-v2
    % 4: delta- >= v2-v1
    % 5: c'v1 = f1 (deletion strain 1)
    % 6: c'v2 = f2 (deletion strain 2)
    
    A = [modelDel1.S sparse(nMets1,nRxns2+2*nCommon);
         sparse(nMets2,nRxns1) modelDel2.S sparse(nMets2,2*nCommon);
         createDeltaMatchMatrix(modelDel1.rxns,modelDel2.rxns);
         modelDel1.c' sparse(1,nRxns2+2*nCommon);
         sparse(1,nRxns1) modelDel2.c' sparse(1,2*nCommon)];

    % Construct the RHS vector
    b = [zeros(nMets1+nMets2+2*nCommon,1);objValDel1 ;objValDel2];

    % Construct the objective (sum of all delta+ and delta-)
    c = [zeros(nRxns1+nRxns2,1);ones(2*nCommon,1)];

    % Construct the ub/lb
    % delta+ and delta- are in [0 10000]
    lb = [modelDel1.lb;modelDel2.lb;zeros(2*nCommon,1)];
    ub = [modelDel1.ub;modelDel2.ub;10000*ones(2*nCommon,1)];

    % Construct the constraint direction vector (G for delta's, E for
    % everything else)
    csense(1:(nMets1+nMets2)) = 'E';
    csense((nMets1+nMets2)+1:(nMets1+nMets2+2*nCommon)) = 'G';
    if (strcmp(osenseStr,'max'))
        csense(end+1) = 'G';
        csense(end+1) = 'G';
    else
        csense(end+1) = 'L';
        csense(end+1) = 'L';
    end
    
    % First find optimal solution value subject to a LP objective 
    if (verbFlag)
        fprintf('STEP 1: Solving linear MOMA upgraded for double knockouts: %d constraints %d variables ',size(A,1),size(A,2));
    end
    
    % Solve the linearMOMA problem
    [LPproblem.A,LPproblem.b,LPproblem.c,LPproblem.lb,LPproblem.ub,LPproblem.csense,LPproblem.osense] = deal(A,b,c,lb,ub,csense,1);
    LPsolution = solveCobraLP(LPproblem);

    if (verbFlag)
        fprintf('%f seconds\n',LPsolution.time);
    end
    
    if (LPsolution.stat > 0)
        solution_full = LPsolution.full
    end
    
    % Use solution_full as an additional constraint for sparseLP
    % Extend the LHS matrix
    % For variable x = [v1;v2;delta+;delta-]
    % Rows:
    % 1: Sdel1*v1 = 0 for the deletion strain 1
    % 2: Sdel2*v2 = 0 for the deletion strain 2
    % 3: delta+ >= v1-v2
    % 4: delta- >= v2-v1
    % 5: c'v1 = f1 (deletion strain 1)
    % 6: c'v2 = f2 (deletion strain 2)
    % 7: c'x = c'*solution_full 
    
    A_0 = [A; c']

    % Extend the RHS vector
    b_0 = [b: c'*solution_full];

    % Construct the ub/lb
    % delta+ and delta- are in [0 10000]
    lb_0 = lb;
    ub_0 = ub;

    % Construct the constraint direction vector (G for delta's, E for
    % everything else)
    csense_0 = [csense; 'E']
    
    %
    if (verbFlag)
        fprintf('STEP 2: Solving linear MOMA upgraded for double knockouts: %d constraints %d variables ',size(A,1),size(A,2));
    end
    
    % Solve the sparseMOMA problem
    [constraint.A,constraint.b,constraint.lb,constraint.ub,constraint.csense] = deal(A_0,b_0,lb_0,ub_0,csense_0);
    
    % Try all non-convex approximations of zero norm and take the best result
    approximations = {'cappedL1','exp','log','SCAD','lp-','lp+'};
    bestResult = n;
    bestAprox = '';
    for i=1:length(approximations)
        solution = sparseLP(char(approximations(i)),constraint);
        if solution.stat == 1
            if bestResult > length(find(abs(solution.x)>eps))
                bestResult=length(find(abs(solution.x)>eps));
                bestAprox = char(approximations(i));
                solutionL0 = solution;
            end
        end
    end
    
    %TO BE CONTINUED
    
    
    
    
    % Call the sparse LP solver
    zeroNormApprox = 'cappedL1'; %This is default apporximation used in optimizeCbModel

    LPsolution = sparseLP(zeroNormApprox,LPproblem);

    if (verbFlag)
%         fprintf('%f seconds\n',LPsolution.time);
%         fprintf('%d number of iterations \n',nIterations);
        fprintf('%Best approximation : %f \n',solution.bestAprox);
    end

    if (LPsolution.stat > 0)
        solutionDel1.x = LPsolution.x(1:nRxns1);
        solutionDel1.f = sum(modelDel1.c.*solutionDel1.x);
        solutionDel2.x = LPsolution.x((nRxns1+1):(nRxns1+nRxns2));
        solutionDel2.f = sum(modelDel2.c.*solutionDel2.x);
%         totalFluxDiff = LPsolution.x;
    end

    if (LPsolution.stat > 0 && minFluxFlag)
        A = [modelDel1.S sparse(nMets1,nRxns2+2*nCommon+2*nRxns1+2*nRxns2);
            sparse(nMets2,nRxns1) modelDel2.S sparse(nMets2,2*nCommon+2*nRxns1+2*nRxns2);
            createDeltaMatchMatrix(modelDel1.rxns,modelDel2.rxns) sparse(2*nCommon,2*nRxns1+2*nRxns2);
            speye(nRxns1,nRxns1) sparse(nRxns1,nRxns2) sparse(nRxns1,2*nCommon) speye(nRxns1,nRxns1) sparse(nRxns1,nRxns1+2*nRxns2);
            -speye(nRxns1,nRxns1) sparse(nRxns1,nRxns2) sparse(nRxns1,2*nCommon) sparse(nRxns1,nRxns1) speye(nRxns1,nRxns1) speye(nRxns1,2*nRxns2);
            sparse(nRxns2,nRxns1) speye(nRxns2,nRxns2) sparse(nRxns2,2*nCommon) sparse(nRxns2,2*nRxns1) speye(nRxns2,nRxns2) sparse(nRxns2,nRxns2);
            sparse(nRxns2,nRxns1) -speye(nRxns2,nRxns2) sparse(nRxns2,2*nCommon) sparse(nRxns2,2*nRxns1) sparse(nRxns2,nRxns2) speye(nRxns2,nRxns2);
            modelDel1.c' sparse(1,nRxns2+2*nCommon+2*nRxns1+2*nRxns2);
            sparse(1,nRxns1) modelDel2.c' sparse(1,2*nCommon+2*nRxns1+2*nRxns2);
            sparse(1,nRxns1+nRxns2) ones(1,2*nCommon) sparse(1,2*nRxns1+2*nRxns2)];
        % Construct the RHS vector
        b = [zeros(nMets1+nMets2+2*nCommon+2*nRxns1+2*nRxns2,1);objValDel1;objValDel2;ceil(totalFluxDiff/tol)*tol];

        % Construct the objective (sum of all delta+ and delta-)
        c = [zeros(nRxns1+nRxns2+2*nCommon,1);ones(2*nRxns1+2*nRxns2,1)];

        % Construct the ub/lb
        % delta+ and delta- are in [0 10000]
        lb = [modelDel1.lb;modelDel2.lb;zeros(2*nCommon+2*nRxns1+2*nRxns2,1)];
        ub = [modelDel1.ub;modelDel2.ub;10000*ones(2*nCommon+2*nRxns1+2*nRxns2,1)];
        csense(1:(nMets1+nMets2)) = 'E';
        csense((nMets1+nMets2)+1:(nMets1+nMets2+2*nCommon+2*nRxns1+2*nRxns2)) = 'G';
        if (strcmp(osenseStr,'max'))
            csense(end+1) = 'G';
            csense(end+1) = 'G';
        else
            csense(end+1) = 'L';
            csense(end+1) = 'L';
        end
        csense(end+1) = 'L';

        if (verbFlag)
            fprintf('Minimizing MOMA flux distribution norms: %d constraints %d variables ',size(A,1),size(A,2));
        end
        csense
        [LPproblem.A,LPproblem.b,LPproblem.lb,LPproblem.ub,LPproblem.csense] = deal(A,b,lb,ub, csense);
        
        [LPsolution, nIterations, bestApprox] = sparseLP(zeroNormApprox,LPproblem);

        if (verbFlag)
            fprintf('%f seconds\n',LPsolution.time);
            fprintf('%d number of iterations \n',nIterations);
            fprintf('%Best approximation : %f \n',bestApprox);
        end

        if (LPsolution.stat > 0)
            solutionDel1.x = LPsolution.x(1:nRxns1);
            solutionDel1.f = sum(modelDel1.c.*solutionDel1.x);
            solutionDel2.x = LPsolution.x((nRxns1+1):(nRxns1+nRxns2));
            solutionDel2.f = sum(modelDel2.c.*solutionDel2.x);
%             totalFluxDiff = LPsolution.obj;
        end

    end
    solutionDel1.stat = LPsolution.stat;
    solutionDel2.stat = LPsolution.stat;
    solStatus = LPsolution.stat;

elseif solutionDel1.stat <= 0 
    warning('Deletion 1 strain FBA problem is infeasible or unconstrained');
    solStatus = solutionDel1.stat;
%     totalFluxDiff = []; 
else 
    warning('Deletion 2 strain FBA problem is infeasible or unconstrained');
    solStatus = solutionDel2.stat;
%     totalFluxDiff = []; 
end


