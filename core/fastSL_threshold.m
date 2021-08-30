function fastSL_threshold(model, threshold, cutoff, order, eliList, atpm)                                  
%% fastSL(model,cutoff,order,eliList,atpm)
% Requires the openCOBRA toolbox
% http://opencobra.sourceforge.net/openCOBRA/Welcome.html
% 
% INPUT
% model (the following fields are required - others can be supplied)       
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         Reaction Names
%OPTIONAL
% threshold      threshold flux for active/non-zero reaction. Default is 0.
% cutoff         cutoff percentage value for lethality.Default is 0.01.
% order          Order of SLs required.Default order is 2. Max value 3.
% eliList        List of reactions to be ignored for lethality
%                analysis:Exchange Reactions, ATPM etc.
% atpm           ATPM Reaction Id in model.rxns if other than 'ATPM'
% 
%
% Aditya Pratapa        6/26/14. 
% N Sowmya Manojna      3 July, 2021

%%
%initCobraToolbox
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

if exist('order', 'var')
    if isempty(order)
        order = 2;
    else
        if (order>3)
        err = MException('ResultChk:OutOfRange', ...
        'Resulting value is outside expected range. Maximum value is 3.');
         throw(err)
        end
    end
else
    order = 2;
end


%Please change this according to your model
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

fname=strcat(model.description,'_Rxn_lethals.mat');


%%

switch order
    case 1
        [Jsl] = singleSL_threshold(model, threshold, cutoff, eliList, atpm);
       
        fprintf('\nSaving Single Lethal Reactions List... ');
        save(fname,'Jsl');
        fprintf('Done. \n');
    case 2
        [Jsl,Jdl] = doubleSL_threshold(model, threshold, cutoff, eliList, atpm);
     
        fprintf('\nSaving Single and Double Lethal Reactions List... ');
        save(fname,'Jsl');
        save(fname,'Jdl','-append');
        fprintf('Done. \n');
    case 3
        [Jsl,Jdl,Jtl]=tripleSL_threshold(model, threshold, cutoff, eliList, atpm);
      
        fprintf('\nSaving Single, Double and Triple Lethal Reactions List... ');
        save(fname,'Jsl');
        save(fname,'Jdl','-append');
        save(fname,'Jtl','-append');
        fprintf('Done. \n');
end

