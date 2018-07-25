function get_output_from_castle(Castle, output_dir)
%GET_CSV_FROM_CASTLE generates output files for mat object Castle
%   Detailed explanation goes here
%   Input:
%   Castle : The structure with models and synthetic lethals data    
%   output_dir : path to output folder to store all the csv files

% Generate Synthetic Lethal xls files in output directory
for iModel = 1:length(Castle.data)
   model_name = Castle.data(iModel).model_name;
   output_dir_i = horzcat(output_dir, '\' ,model_name);
   
   Jsl = Castle.data(iModel).Jsl;
   Jdl = Castle.data(iModel).Jdl;
   
   data1 = {'Single Lethal Reactions', Jsl{:}};
   data2 = {'Double Lethal Reactions','Rxn A', Jdl{:,1};
            '','Rxn B', Jdl{:,2}};
end    
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



end

