function [ Rerouting_List ] = Rerouting_RxnList( model, rxnList , FVA_cutoff)
%[Rerouting_List] = Rerouting_rxnList(model, rxnList , FVA_cutoff)
%   Finds reactions in minimal rerouting when given set of
%reactions are deleted
% INPUT
% model - SBML model
% rxnList - List of reactions
% Optional INPUT
% FVA_cutoff - Forced lower bound for the reactions (WT solution) Default:0.5 of FVA max
% OUTPUT
% Rerouting_List - Structure rerturning set of reactions in minimal rerouting
% Omkar 21 Oct, 2017

if exist('FVA_cutoff', 'var')
    if isempty(cutoff)
        FVA_cutoff = 0.5;
    end
else
    FVA_cutoff = 0.5;
end
cutOff=1e-4; % Rerouting

[minFlux,maxFlux] = fluxVariability(model,0,'max',rxnList);
rxnList_lb=FVA_cutoff*maxFlux;

rxnListID=find(ismember(model.rxns,rxnList));
minRerouting(length(rxnList)).rxns=[];
minRerouting(length(rxnList)).diff=[];

for iRxn=2532:length(rxnList)
    % Find WT flux distribution so that iRxn has minimal flux of rxnList_lb
    if rxnList_lb(iRxn)>0
        model_i=model;
        model_i.lb(rxnListID(iRxn))=rxnList_lb(iRxn);
        solWT_i=optimizeCbModel(model_i,'max','zero');
        fluxWT_i=solWT_i.x;
        
        % Delete the reaction to find rerouting
        model_i.lb(iRxn)=0; model_i.ub(iRxn)=0;
        solKO_i=fluxMOMAzn(model_i,fluxWT_i);
        fluxMT=solKO_i.x;
        diff=fluxMT-fluxWT_i;
        
        minRerouting(iRxn).rxns=model_i.rxns(find(abs(diff)>cutOff));
        minRerouting(iRxn).diff=diff(find(abs(diff)>cutOff));
    end    
end

fname=sprintf('minReRxn.csv');
     fid = fopen(fullfile(fname),'wt');
     if fid>0
         for k=1:length(minRerouting)
             fprintf(fid,'%s,%d,',char(model.rxns(k)),length(minRerouting(k).rxns));
             for m=1:length(minRerouting(k).rxns)
                 fprintf(fid,'%s,',char(minRerouting(k).rxns(m)));
             end
             fprintf(fid,'\n , %d,',length(minRerouting(k).rxns));   
             for m=1:length(minRerouting(k).rxns)
                 fprintf(fid,'%f,',minRerouting(k).diff(m));
             end
             fprintf(fid,'\n'); 
         end
         fclose(fid);
     end

Rerouting_List=minRerouting;     
end

