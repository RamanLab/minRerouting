% Degree of rerouting
% Generating input for Cytoscape for 6 lethals of different size of
% rerouting. For each lethal pair run the function GetNetwork 
% 1. {'IPDDI','IPDPS'} - 51
% 2. {'G5SD','NACODA'} - 45
% 3. {'3OAS160','FACOAE160'} - 21
% 4. {'TKT2','RPE'} - 127
% 5. {'PGK','O2tex'} - 106
% 6. {'PPKr','PPA'} - 74
iLeth=1;
minSet_i=minRerouted(iLeth).rxns(abs(minRerouteddiff(iLeth,:))>0.0001);
net=getNetwork(minSet_i,model);
 fid = fopen('FVA_DL1.csv','wt');
 fprintf(fid,'source,target,Column 3,Column 4\n');
 if fid>0
     for k=1:size(net,1)
         fprintf(fid,'%s,%s,%s,%s\n',net{k,:});
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
%     if length(minRerouted(IdxLeth(k)).rxns)<50
%         minSetClass_all=[minSetClass_all; minRerouted(IdxLeth(k)).rxns];
%     else 
        minSetClass_all=[minSetClass_all; minRerouted(IdxLeth(k)).rxns(abs(minRerouteddiff(IdxLeth(k),:))>0.0001)];
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
