function get_output_from_castle(Castle, output_dir, norm)
%get_output_from_castle(Castle, output_dir)
% get_output_from_castle generates output files for mat object Castle
%   Input:
%   Castle : The structure with models and synthetic lethals data
%   output_dir : path to output folder to store all the csv files

% Generate Synthetic Lethal xls files in output directory
for iModel = 1:length(Castle.data)
    model_name = Castle.data(iModel).model_name;
    output_dir_i = strcat(output_dir, '/' ,model_name);
    
    Jsl = Castle.data(iModel).Jsl;
    Jdl = Castle.data(iModel).Jdl;
    
    data1 = {'Single Lethal Reactions', Jsl{:}};
    data2 = {'Double Lethal Reactions','Rxn A', Jdl{:,1};
        '','Rxn B', Jdl{:,2}};
    
    fname1 = strcat(output_dir_i, '/', norm, '_SyntheticLethal_list.xlsx');
    %xlswrite(fname1, transpose(data1), 1);
    %xlswrite(fname1, transpose(data2), 2);
    dat1 = cell2table(transpose(data1));
    writetable((dat1),fname1,'sheet',1);
    dat2 = cell2table(transpose(data2));
    writetable(dat2, fname1, "sheet", 2);
    %writecell(transpose(data1),"SL_sl.csv")
    %writecell(transpose(data2),"SL_dl.csv")

    % MinReSet Excel file
    minReSets = Castle.data(iModel).minRe;
    data3 = {};
    data3(1,:) = {'Synthetic Lethal Pair', '', 'Size of minReSet','Total flux difference(abs)','short', 'long', 'common', 'Reactions in minReSet/ Flux difference'};
    for k = 1:length(Jdl)
        data3(2*k,1:7) = {Jdl{k,1},Jdl{k,2},length(minReSets(k).rxns), sum(abs(minReSets(k).diff_flux)), length(minReSets(k).PathShort), length(minReSets(k).PathLong), length(minReSets(k).pathCommon)};
        
        for m=1:length(minReSets(k).rxns)
            data3(2*k,7+m) = minReSets(k).rxns(m);
            data3(2*k+1,7+m) = num2cell(minReSets(k).diff_flux(m));
        end
    end
    fname3 = strcat(output_dir_i, '/', norm, '_minReroutingSets.xlsx');
    %xlswrite(fname3, data3, 1);
    dat3 = cell2table(data3);
    writetable(dat3, fname3)
    %writecell(data3, "TRIAL_one_minrerouting.csv")
end
end

