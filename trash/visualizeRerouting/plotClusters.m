function [ Cyto ] = plotClusters( model,Cluster, Jdl, minRerouting )
%PLOTCLUSTERS generates CSV files for all clusters identified to be used in
%cytoscape
for i=1:length(Cluster)
    Leth_cat=strcat(Cluster(i).Jdl(:,1),Cluster(i).Jdl(:,2));
    minSetClass_all=[];

    Jdl_cat=strcat(Jdl(:,1),Jdl(:,2));
    IdxLeth=[];
    for m=1:length(Leth_cat)
        IdxLeth=[IdxLeth; find(strcmp(Leth_cat(m),Jdl_cat))];
    end    
    for k=1:length(IdxLeth)
            minSetClass_all=[minSetClass_all; minRerouting(IdxLeth(k)).rxns];
    end

    minSet_i=unique(minSetClass_all);
    net=getNetwork(minSet_i,model);
    fname=sprintf('clust_%d.csv',i)
     fid = fopen(fullfile(fname),'wt');
     fprintf(fid,'source,target,Column 3\n');
     if fid>0
         for k=1:size(net,1)
             fprintf(fid,'%s,%s,%s\n',net{k,:});
         end
         fclose(fid);
     end
     Cyto(i).net=net;
end
end

