function [Cluster,nClust] = clusterMinSet( model, Jdl, minRerouted )
%Clustering of lethals based on minRerouting sets
p=zeros(length(Jdl),length(model.rxns));
for iLeth=1:length(Jdl)
    if ~isempty(minRerouted(iLeth).rxns)
        p(iLeth,find(ismember(model.rxns,minRerouted(iLeth).rxns)))=1;
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
    fprintf('Cluster %d:\n',k)
    Cluster(m).Jdl=Jdl(c==k,:); m=m+1;
    Jdl(c==k,:)
    unique(sum(p(find(c==k),:),2))
end
nClust=length(Cluster);
end

