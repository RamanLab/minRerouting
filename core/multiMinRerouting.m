function [Castle] = multiMinRerouting(Castle)
%MULTIMINREROUTING Summary of this function goes here
%   Detailed explanation goes here
nModels = length(Castle.data);
for j= 1:nModels
    disp("Finding set of minimal rerouting reactions for double lethals in " + Castle.data(j).model_name);
%     Castle.data(j).minRe =
%     minReroutingRxns(Castle.data(j).model,Castle.data(j).Jdl); % Old
%     version
%    Castle.data(j).minRe = minReroutingRxns_l2_moma(Castle.data(j).model,Castle.data(j).Jdl);    %quadraticMOMA version
%    Castle.data(j).minRe = minReroutingRxns_l1_moma(Castle.data(j).model,Castle.data(j).Jdl);    %linearMOMA version
   Castle.data(j).minRe = minReroutingRxns_l0_moma(Castle.data(j).model,Castle.data(j).Jdl);    %sparseMOMA version

end    
end

