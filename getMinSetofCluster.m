function [ minReSet ] = getMinSetofCluster(lethSet,Jdl,minRerouted)
%GETMINSETOFCLUSTER Summary of this function goes here
%   Detailed explanation goes here

    Leth_cat=strcat(lethSet(:,1),lethSet(:,2));
    minSetClass_all=[];

    Jdl_cat=strcat(Jdl(:,1),Jdl(:,2));
    IdxLeth=[];
    for m=1:length(Leth_cat)
        IdxLeth=[IdxLeth; find(strcmp(Leth_cat(m),Jdl_cat))];
    end    
    for k=1:length(IdxLeth)
            minSetClass_all=[minSetClass_all; minRerouted(IdxLeth(k)).rxns];
    end
    
    minReSet=unique(minSetClass_all);
end

