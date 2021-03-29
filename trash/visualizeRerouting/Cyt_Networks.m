% Degree of rerouting
% Generating input for Cytoscape for 6 lethals of different size of
% rerouting. For each lethal pair run the function GetNetwork 
% 1. {'NI2tpp','NI2uabcpp'} - 212
% 2. {'G3PAT160','APG3PAT160'} -138
% 3. {'G5SD','NACODA'} - 141
% 4. {'I4FE4ST','S4FE4ST'} - 204
% 5. {'AACPS7','3HAD120'} - 7
% 6. {'PGM','O2tex'} - 250 
iLeth=250;
minSet_i=minReFBA(iLeth).cluster;
net=getNetwork(minSet_i,model);
fname=sprintf('DLPairFBA_%d.csv',iLeth)
     fid = fopen(fullfile(fname),'wt');
     fprintf(fid,'source,target,Column 3\n');
     if fid>0
         for k=1:size(net,1)
             fprintf(fid,'%s,%s,%s\n',net{k,:});
         end
         fclose(fid);
     end
% Classification of lethals
% 1. Acetyl pathways - 21 Lethals
Leth_cat=strcat(LethalsAcetylPathway(:,1),LethalsAcetylPathway(:,2));

% 2. NDPK Pathway - 16 Lethals
Leth_cat=strcat(LethalsNDPKPathway(:,1),LethalsNDPKPathway(:,2));

% 3. TKT Pathway - 12 lethals
Leth_cat=strcat(LethalsTKT(:,1),LethalsTKT(:,2));

% 4. Fe Transport Pathway - 10
Leth_cat=strcat(LethalsFeTrans(:,1),LethalsFeTrans(:,2));

% 5. Diphosphate Pathway - 3 
Leth_cat=strcat(LethalsDiPh(:,1),LethalsDiPh(:,2));

% 6. NAD Transport - 6
Leth_cat=strcat(LethalsNAD(:,1),LethalsNAD(:,2));

% 7. Others - 15
Leth_cat=strcat(LethalsOthers(:,1),LethalsOthers(:,2));

% 7. O2 Transport(a) - 7
Leth_cat=strcat(LethalsO2Tran(:,1),LethalsO2Tran(:,2));


minSetClass_all=[];

Jdl_cat=strcat(Jdl(:,1),Jdl(:,2));
IdxLeth=[];
for m=1:length(Leth_cat)
    IdxLeth=[IdxLeth; find(strcmp(Leth_cat(m),Jdl_cat))];
end    
for k=1:length(IdxLeth)
%     if length(minRerouting(IdxLeth(k)).rxns)<50
%         minSetClass_all=[minSetClass_all; minRerouting(IdxLeth(k)).rxns];
%     else 
        minSetClass_all=[minSetClass_all; minRerouting(IdxLeth(k)).rxns(abs(minReroutingdiff(IdxLeth(k),:))>0.0001)];
%     end
end

minSet_i=unique(minSetClass_all);
net=getNetwork(minSet_i,model);
 fid = fopen('pathOthers.csv','wt');
 fprintf(fid,'source,target,Column 3,Column 4\n');
 if fid>0
     for k=1:size(net,1)
         fprintf(fid,'%s,%s,%s,%s\n',net{k,:});
     end
     fclose(fid);
 end
