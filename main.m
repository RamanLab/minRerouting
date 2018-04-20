function minReOutputs = minRe(filename, SL_data)
%% minReOutputs = minRe(filename)
% Input :
% filename - path of xml model file 
% Optional
% SL_data  - FastSL output file (With structure including at least Jsl and Jdl) ; if empty run fastSL

initCobraToolbox;
addpath('core\findRerouting','core\visualizeRerouting'); %adding libraries
delete('Output\*'); % delete output folder to save newly genrated data for current model

model = readCbModel(filename);

if (nargin <2 || isempty(cutOff))
    eliList=model.rxns(find(findExcRxns(model)));
    fastSL(model,0.05,2,eliList);
    SL_dir = dir('*_Rxn_lethals.mat');
    load(SL_dir.name)
    movefile (SL_dir.name,'Output')
else
    load(SL_data);
    movefile(SL_data,'Output')
end

% Generate Synthetic Lethal xls files in output directory
fname = 'Output\SyntheticLethal_list.xlsx';
data1 = {'Single Lethal Reactions', Jsl{:}};
xlswrite(fname, transpose(data1), 1);
data2 = {'Double Lethal Reactions','Rxn A', Jdl{:,1};
            '','Rxn B', Jdl{:,2}};
xlswrite(fname, transpose(data2), 2);

% Generate min rerouting data in output folder
minReSets=minReroutingRxns(model,Jdl,0.0001,3,'True');

for k=1:length(minReSets)
    size_minReSets(k,1) = length(minReSets(k).rxns);
    size_minReShort(k,1) = length(minReSets(k).PathShort);
    size_minReLong(k,1) = length(minReSets(k).PathLong);
end    
% Histogram plot
hist(size_minReSets,15)
title('Distribution of degree of rerouting across double lethals','fontSize',20)
xlabel('Number of reactions being rerouted','fontSize',20)
ylabel('Number of double lethals pairs','fontSize',20)
print('Output\hist_minReSets','-dpng')

% Scatter plot
[size_uniq]=unique([size_minReShort, size_minReLong],'rows');
tbl=tabulate([1000*size_minReShort+size_minReLong]);

b1=size_uniq(:,1)\size_uniq(:,2);
yCalc1 = b1*size_uniq(:,1);
figure(3)
scatter(size_uniq(:,1),size_uniq(:,2),20*tbl(:,2))
hold on
plot(size_uniq(:,1),yCalc1)
title('Comparison between lengths of alternative paths','fontSize',20)
xlabel('Number of reactions in Shorter Path','fontSize',20)
ylabel('Number of reactions in Longer Path','fontSize',20)
print('Output\scatter_pathlen','-dpng')




    fprintf(fid,'%s,%s,%d,',Jdl{k,1},Jdl{k,2},length(minRerouting(k).rxns));
    for m=1:length(minRerouting(k).rxns)
        fprintf(fid,'%s,',minRerouting(k).rxns{m});
    end
    fprintf(fid,'\n');
end

fname=sprintf('minReAll_iML.csv');
     fid = fopen(fullfile(fname),'wt');
     if fid>0
         for k=1:length(minRerouting)
             fprintf(fid,'%s,%s,%d,',Jdl{k,1},Jdl{k,2},length(minRerouting(k).rxns));
             for m=1:length(minRerouting(k).rxns)
                 fprintf(fid,'%s,',minRerouting(k).rxns{m});
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
             fprintf(fid,'%d,%d,',k,size(Clusters(k).Jdl,1));
             
             for m=1:size(Clusters(k).Jdl,1)
                 fprintf(fid,'%s+%s,',char(Clusters(k).Jdl(m,1)),char(Clusters(k).Jdl(m,2)));
             end
             fprintf(fid,'\n');         
         end
         fclose(fid);
     end
     

lethSet_1=ClusterEcoli(1).Jdl;

Cyto=plotClusters(model,Clusters,Jdl,minRerouting);

end