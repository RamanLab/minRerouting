function [ Classified ] = classifyLethals(model, Jdl, sol_x)
%UNTITLED3 Summary of this function goes here
%   sol_x  FBA solution
Classified.Redundant=[];
Classified.Plastic=[];

for i=1:length(Jdl)
    if abs(sol_x(find(ismember(model.rxns,Jdl(i,1)))))>10e-6 && abs(sol_x(find(ismember(model.rxns,Jdl(i,2)))))>10e-6
        Classified.Redundant=[Classified.Redundant; Jdl(i,:)];
    else 
        Classified.Plastic=[Classified.Plastic; Jdl(i,:)];
    end
end
end