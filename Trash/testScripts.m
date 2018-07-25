% minRerouting
%Package to understand the compensatory mechanisms performing similar metabolic functions through analysis of Synthetic Lethals

% Alternate Routes in Networks
%The rerouting of the reaction fluxes in the metabolic network to produce an essential function seems to be non-trivial and involves complex structrues. This package helps to identify minimal rerouting corresponding synthetic lethal pairs using constraint based tools like MOMA.   

% Steps in analysis
%1. Get FastSL toolbox from https://github.com/RamanLab/FastSL

%Here we explain the example of some of the models of pathogenic organims downloaded from BIGG database. Exchange reactions are exluded from the analysis to focus on intracellular properties of metabolic networks

%2. Run FastSL on COBRA model/s object to identify synthetic lethal pairs for given models.
% Output structure Castle represents Comparitice Analysis of SynThetic LEthals
addpath('core\findRerouting\','core\compareLethals\','core\visualizeRerouting')
model_names = {'iIT341', 'iML1515', 'iNJ661', 'iPC815', 'iYL1228', 'STM_v1_0'} % Recommneded to rename model.mat file to modelNames.mat file
path_to_models = {'examples/iIT341/', 'examples/iML1515/', 'examples/iNJ661/','examples/iPC815/', 'examples/iYL1228/', 'examples/STM_v1_0/'}
Castle_tmp2 = getFastSL(model_names, path_to_models, 2); % 2 for double lethals and 3 for triple lethals
% To Do: modify for without fastSL run in cases where SL are identified 

%3. Comparitive analysis synthetic lethals 
% Use multiData input Struct to run compareLethals
Castle = compareLethals(Castle, 2); 
% To Do: modify for TL compare





















% Test Scripts
% This script is used to test and analyse performance of minRe algorithm
% Trying to resolve any inconsistencies observed in algorithm is main
% function of this script

% Implementing relative tolerance
% example SL:  G3PD2, GLYK (83) 
iLeth = 83;
delIdx_1=find(strcmp(Jdl(iLeth,1),model.rxns)); delIdx_2=find(strcmp(Jdl(iLeth,2),model.rxns));
modeldel_1=model; modeldel_2=model;    
modeldel_1.lb(delIdx_1)=0;    modeldel_1.ub(delIdx_1)=0;
modeldel_2.lb(delIdx_2)=0;    modeldel_2.ub(delIdx_2)=0;
    Iter=6;
V1=zeros(length(model.rxns),Iter); V2=zeros(length(model.rxns),Iter);

sol_11=optimizeCbModel(modeldel_1,'max','zero'); %Sol_11 is MT1 Iteration 1
V1(:,1)=sol_11.x;
fluxDel = sol_11.x;
    
    sol_21=Copy_of_fluxMOMAzn(modeldel_2,V1(:,1));   %Sol_21 is MT2 Iteration 1
    V2(:,1)=sol_21.x;
    V2_stat(1,1) = sol_21.stat;



% 1. Number of iterations to be run before calculating minRerouting sets
% Current algorithm is ran 3 iterations by defaault. This may not be
% optimal and needs to be tested 
% minRerouting_5=minReroutingRxns(model,Jdl,0.0001,5); % Uncomment the inconsistency 1 part of core script of minRe
for k = 1:length(minRerouting_5_delta01)
    size_iter_data(k,1) = length(minRerouting_5_delta01(k).Diff1);
    size_iter_data(k,2) = length(minRerouting_5_delta01(k).Diff2);
    size_iter_data(k,3) = length(minRerouting_5_delta01(k).Diff3);
    size_iter_data(k,4) = length(minRerouting_5_delta01(k).Diff4);
    size_iter_data(k,5) = length(minRerouting_5_delta01(k).Diff5);
end    
plot(size_iter_data(:,3:5))
xlswrite('Iterations.xlsx',size_iter_data)

% Exploring the incosistent lethal pairs
% G3PD2, GLYK
modeldel_1=model;
modeldel_2=model;
iLeth = 83;

delIdx_1=find(strcmp(Jdl(iLeth,1),model.rxns));
delIdx_2=find(strcmp(Jdl(iLeth,2),model.rxns));
    
    modeldel_1.lb(delIdx_1)=0;
    modeldel_1.ub(delIdx_1)=0;
    
    modeldel_2.lb(delIdx_2)=0;
    modeldel_2.ub(delIdx_2)=0;
    cutOff = 0.000001;  delta = 0.2; Iter = 6;
    V1=zeros(length(model.rxns),Iter);
    V2=zeros(length(model.rxns),Iter);

    sol_11=optimizeCbModel(modeldel_1,'max','zero'); %Sol_11 is MT1 Iteration 1
    sol_11.x(sol_11.x<cutOff)=0;
    V1(:,1)=sol_11.x;
    
    V1_stat(1,1)=sol_11.stat;
    sol_21=fluxMOMAzn(modeldel_2,V1(:,1));   %Sol_21 is MT2 Iteration 1
    sol_21.x(sol_21.x<cutOff)=0
    V2(:,1)=sol_21.x;
    
    for j=2:Iter
        sol_1j=fluxMOMAzn(modeldel_1,V2(:,j-1));
        sol_1j.x(sol_1j.x<cutOff)=0
        V1(:,j)=sol_1j.x;
        
        sol_2j=fluxMOMAzn(modeldel_2,V1(:,j));        
        sol_2j.x(sol_2j.x<cutOff)=0
        V2(:,j)=sol_2j.x;
    end

    Diff_mat = abs(V1-V2);
    diff.Diff1 = find(Diff_mat(:,1)>delta*V1(:,1) & Diff_mat(:,1)>cutOff);
    diff.Diff2 = find(Diff_mat(:,2)>delta*V1(:,2) & Diff_mat(:,2)>cutOff);
    diff.Diff3 = find(Diff_mat(:,3)>delta*V1(:,3) & Diff_mat(:,3)>cutOff);
    diff.Diff4 = find(Diff_mat(:,4)>delta*V1(:,4) & Diff_mat(:,4)>cutOff);
    diff.Diff5 = find(Diff_mat(:,5)>delta*V1(:,5) & Diff_mat(:,5)>cutOff);
    diff.Diff6 = find(Diff_mat(:,6)>delta*V1(:,6) & Diff_mat(:,6)>cutOff);

    
 net1=getNetwork(diff_rxn.Diff1,model);
 net2=getNetwork(diff_rxn.Diff2,model);
 
 fname=sprintf('net_2.csv')
     fid = fopen(fullfile(fname),'wt');
     fprintf(fid,'source,target,Column 3\n');
     if fid>0
         for k=1:size(net2,1)
             fprintf(fid,'%s,%s,%s\n',net2{k,:});
         end
         fclose(fid);
     end
   
 Cyto=plotClusters(model,Clusters,Jdl,minRerouting);

    