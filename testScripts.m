% Test Scripts
% This script is used to test and analyse performance of minRe algorithm
% Trying to resolve any inconsistencies observed in algorithm is main
% function of this script

% 1. Number of iterations to be run before calculating minRerouting sets
% Current algorithm is ran 3 iterations by defaault. This may not be
% optimal and needs to be tested 
% minRerouting_5=minReroutingRxns(model,Jdl,0.0001,5); % Uncomment the inconsistency 1 part of core script of minRe
for k = 1:length(minRerouting_5)
    size_iter_data(k,1) = length(minRerouting_5(k).Diff1);
    size_iter_data(k,2) = length(minRerouting_5(k).Diff2);
    size_iter_data(k,3) = length(minRerouting_5(k).Diff3);
    size_iter_data(k,4) = length(minRerouting_5(k).Diff4);
    size_iter_data(k,5) = length(minRerouting_5(k).Diff5);
end    
plot(size_iter_data(:,3:5))
xlswrite('Iterations.xlsx',size_iter_data)