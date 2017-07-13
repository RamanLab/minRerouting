initCobraToolbox;
% STM_v1.0
load('Sty_STM_v1.0.mat')
load('STM_v1.0_Rxn_lethals.mat')

minRerouting=minReroutingRxns(model,Jdl,0.0001);

fname=sprintf('minReRxns_STM.csv');
     fid = fopen(fullfile(fname),'wt');
     if fid>0
         for k=1:length(Jdl)
             fprintf(fid,'%s,%s,%d,',char(Jdl(k,1)),char(Jdl(k,2)),length(minRerouting(k).rxns));
             for m=1:length(minRerouting(k).rxns)
                 fprintf(fid,'%s,',char(minRerouting(k).rxns(m)));
             end
             fprintf(fid,'\n');         
         end
         fclose(fid);
     end
     
Clusters=clusterMinSet(model,Jdl,minRerouting);
fname=sprintf('Clusters_Auto.csv');
     fid = fopen(fullfile(fname),'wt');
     if fid>0
         for k=1:length(Clusters)
             fprintf(fid,'%d,%d,',k,length(Clusters(k).Jdl));
             for m=1:length(Clusters(k).Jdl)
                 fprintf(fid,'%s+%s,',char(Clusters(k).Jdl(m,1)),char(Clusters(k).Jdl(m,2)));
             end
             fprintf(fid,'\n');         
         end
         fclose(fid);
     end

Cyto=plotClusters(model,Clusters,Jdl,minRerouting);
