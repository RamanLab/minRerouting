function [Jsl] = singleSL_threshold(model, threshold, cutoff, eliList, atpm)
%% [Jsl] = singleSL(model,cutoff,eliList,atpm)
% INPUT
% model (the following fields are required - others can be supplied)       
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         Reaction Names
% OPTIONAL
% threshold      threshold flux for active/non-zero reaction. Default is 0.
% cutoff         cutoff percentage value for lethality.Default is 0.01.
% eliList        List of reactions to be ignored for lethality
%                analysis:Exchange Reactions, ATPM etc.
% atpm           ATPM Reaction Id in model.rxns if other than 'ATPM'
% OUTPUT
% Jsl            Single lethal reactions identified
% 
% Aditya Pratapa        6/26/14. 
% N Sowmya Manojna      3 July, 2021

if exist('threshold', 'var')
    if isempty(threshold)
        threshold = 0;
    end
else
    threshold = 0;
end

if exist('cutoff', 'var')
    if isempty(cutoff)
        cutoff = 0.01;
    end
else
    cutoff = 0.01;
end

if exist('atpm', 'var')
    if isempty(atpm)
        atpm = 'ATPM'; %Reaction Id of ATP maintenance reaction- by default it takes 'ATPM'
    end
else
    atpm = 'ATPM';
end

if exist('eliList', 'var')
    if isempty(eliList)
        eliList = model.rxns(ismember(model.rxns,atpm)); %To eliminate ATPM.
    end
else
    eliList = model.rxns(ismember(model.rxns,atpm));
end



Jsl = [];

%Step1 Identify Single Lethal Reactions...
%Identify minNorm flux distribution
solWT = optimizeCbModel(model, 'max', 'one');
grWT = solWT.f;

% Consider greater than threshold than  neq 0
Jnz = find(abs(solWT.x) > threshold);

if (~isempty(eliList))
    % Index of reactions not considered for lethality analysis
    eliIdx = find(ismember(model.rxns, eliList));
    % Jnz
    Jnz = Jnz(~ismember(Jnz, eliIdx));
end
h = waitbar(0,'0.00','Name','Identifying Jsl...');

% Identify Single Lethal Reaction Deletions...
modeldel = model;

for iRxn = 1:length(Jnz)
    delIdx_i = Jnz(iRxn);
    modeldel.lb(delIdx_i) = 0;
    modeldel.ub(delIdx_i) = 0;

    solKO_i = optimizeCbModel(modeldel);
    if (solKO_i.f<cutoff*grWT || isnan(solKO_i.f))
        Jsl = [Jsl;delIdx_i];
    end
   % Reset bounds on idx reaction
    modeldel.lb(delIdx_i) = model.lb(delIdx_i);
    modeldel.ub(delIdx_i) = model.ub(delIdx_i);
    waitbar(iRxn/length(Jnz),h,[num2str(round(iRxn*100/length(Jnz))) '% completed...']);

end
close(h);
Jsl = model.rxns(Jsl);
end