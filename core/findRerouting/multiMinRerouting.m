function [Castle] = multiMinRerouting(Castle)
%MULTIMINREROUTING Summary of this function goes here
%   Detailed explanation goes here
nModels = length(Castle.data);
for j= 1:nModels
    disp("Finding set of minimal rerouting reactions for double lethals in " + Castle.data(j).model_name);
    Castle.data(j).minRe = minReroutingRxns(Castle.data(j).model,Castle.data(j).Jdl); 
end    
end

