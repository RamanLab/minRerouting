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
    
    fname1 = horzcat(output_dir_i, '\SyntheticLethal_list.xlsx');
    xlswrite(fname1, transpose(data1), 1);
    xlswrite(fname1, transpose(data2), 2);
    
    % MinReSet Excel file
    minReSets = Castle.data(iModel).minRe;
    data3 = {};
    data3(1,:) = {'Synthetic Lethal Pair', '', 'Size of minReSet','Total flux difference(abs)', 'Reactions in minReSet/ Flux difference'};
    for k = 1:length(Jdl)
        data3(2*k,1:4) = {Jdl{k,1},Jdl{k,2},length(minReSets(k).rxns), sum(abs(minReSets(k).diff))};
        
        for m=1:length(minReSets(k).rxns)
            data3(2*k,4+m) = minReSets(k).rxns(m);
            data3(2*k+1,4+m) = num2cell(minReSets(k).diff(m));
        end
    end
    fname3 = horzcat(output_dir_i, '\minReroutingSets.xlsx');
    xlswrite(fname3, data3, 1);
end
end

