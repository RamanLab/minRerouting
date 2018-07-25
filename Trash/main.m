function minReOutputs = minRe(filename, SL_data)
%% minReOutputs = minRe(filename)
% Input :
% filename - path of xml model file 
% Optional
% SL_data  - FastSL output file (With structure including at least Jsl and Jdl) ; if empty run fastSL

% initCobraToolbox;
addpath('core\findRerouting','core\visualizeRerouting'); %adding libraries

model = readCbModel(filename);
save('model.mat','model')
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

% Generate min rerouting data in output folder
minReSets=minReroutingRxns(model,Jdl,0.000001,0,3,'True'); 
save('Output\minReSets.mat', 'minReSets'); 

%% Figures 
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

%% Excel files data 
% Generate Synthetic Lethal xls files in output directory
data1 = {'Single Lethal Reactions', Jsl{:}};
data2 = {'Double Lethal Reactions','Rxn A', Jdl{:,1};
            '','Rxn B', Jdl{:,2}};

fname1 = 'Output\SyntheticLethal_list.xlsx';
xlswrite(fname1, transpose(data1), 1);        
xlswrite(fname1, transpose(data2), 2);

% MinReSet Excel file 
data_matrix = {};
data_matrix(1,:) = {'Synthetic Lethal Pair', '', 'Size of minReSet','Total flux difference(abs)', 'Reactions in minReSet/ Flux difference'};
for k = 1:length(Jdl)
    data_matrix(2*k,1:4) = {Jdl{k,1},Jdl{k,2},length(minReSets(k).rxns), sum(abs(minReSets(k).diff))};    

    for m=1:length(minReSets(k).rxns)
        data_matrix(2*k,4+m) = minReSets(k).rxns(m);
        data_matrix(2*k+1,4+m) = num2cell(minReSets(k).diff(m));
    end
end    

fname2= 'Output\minReSetsData.xlsx';
xlswrite(fname2, data_matrix)

% coreMinReSet Excel file 
data_matrix_c = {};
data_matrix_c(1,:) = {'Synthetic Lethal Pair', '', 'Size of coreminReSet','Total flux difference(abs)', 'Reactions in coreminReSet/ Flux difference'};
for k = 1:length(Jdl)
    data_matrix_c(2*k,1:4) = {Jdl{k,1},Jdl{k,2},length(coreReSets(k).rxns), sum(abs(coreReSets(k).diff))};    

    for m=1:length(coreReSets(k).rxns)
        data_matrix_c(2*k,4+m) = coreReSets(k).rxns(m);
        data_matrix_c(2*k+1,4+m) = num2cell(coreReSets(k).diff(m));
    end
end    

fname3= 'Output\coreminReSetsData.xlsx';
xlswrite(fname3, data_matrix_c)

     
Clusters = clusterMinSet(model,Jdl,minReSets);
fname = sprintf('Clusters_Auto.csv');
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