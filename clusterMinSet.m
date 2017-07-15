function [Cluster] = clusterMinSet( model, Jdl, minRerouting )
%Clustering of lethals based on minRerouting sets
p=zeros(length(Jdl),length(model.rxns));
for iLeth=1:length(Jdl)
    if ~isempty(minRerouting(iLeth).rxns)
        p(iLeth,find(ismember(model.rxns,minRerouting(iLeth).rxns)))=1;
    end
end    
D=pdist(p,'jaccard');
l=linkage(D);
c=cluster(l,'Cutoff',eps);
m=1;
for k=1:max(c)
    if(sum(c==k)<=1)
        continue
    end
%     fprintf('Cluster %d:\n',k)
    Cluster_temp(m).Jdl=Jdl(c==k,:); 
    len(m)=length(Jdl(c==k,:)); m=m+1;
end
[Y,I]=sort(len,'descend');
for Idx=1:m-1
    Cluster(Idx).Jdl=Cluster_temp(I(Idx)).Jdl;
end
fprintf('Total number of clusters identified: %d\n',m-1)
end
