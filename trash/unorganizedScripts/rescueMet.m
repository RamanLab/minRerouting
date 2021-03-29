function [RescMet] = rescueMet( model,Jdl,minRerouting )
%rescueMet identifies rescue metabolites based on minReroting sets
% INPUT
% model (the following fields are required - others can be supplied)
%   S            Stoichiometric matrix
%   b            Right hand side = dx/dt
%   c            Objective coefficients
%   lb           Lower bounds
%   ub           Upper bounds
%   rxns         Reaction Names
% Jdl            List of reaction pairs(double lethals generally) for identifying the flux rerouting
% minRerouting   Structure containing details of alternate routes for a DL
%
%OUTPUT
% RescMet        List of rescue metabolites for DLs
%
% Omkar Mohite       13 Jul,2017.

% Idea is to find the common end metabolite in the alternate pathways
for m=1:length(model.mets)
    freq_mets(m,1)=length(find(model.S(m,:)));
end
highFreq_Met=find(freq_mets>10); %freqeuntly occuring metabolite (more than 30 rxns) 
for iLeth=1:length(Jdl)
    Idx_PathShort=find(ismember(model.rxns,minRerouting(iLeth).PathShort));
    Idx_PathLong=find(ismember(model.rxns,minRerouting(iLeth).PathLong));
    met_short=[]; met_long=[];
    for k=1:length(Idx_PathShort)
        met_short=[met_short;find(model.S(:,Idx_PathShort(k))>0)];
    end
    for k=1:length(Idx_PathLong)
        met_long=[met_long;find(model.S(:,Idx_PathLong(k))>0)];
    end
    common_temp=met_long(find(ismember(met_long,met_short)));
    CommonMet(iLeth).mets=common_temp(find(~ismember(common_temp,highFreq_Met))); %Exclude high freq mets
end         
end
