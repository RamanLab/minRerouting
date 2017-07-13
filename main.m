initCobraToolbox;
% load('iAF1260.mat')
% load('iAF1260_Rxn_lethals.mat')
% 
% % BiGG model
% load('iJO1366.mat');
% model=iJO1366;
% eliList=model.rxns(find(findExcRxns(model)));
% fastSL(model,0.01,2,eliList)
% 
% % Model from Paper
% load('iJOPaper.mat');
% model=iJOPaper;
% eliList=model.rxns(find(findExcRxns(model)));
% fastSL(model,0.05,2,eliList)

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
     
Divide=divideMinSet(model,minRerouting,Jdl);     
     
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
     

lethSet_1=ClusterEcoli(1).Jdl;

Cyto=plotClusters(model,Clusters,Jdl,minRerouting);