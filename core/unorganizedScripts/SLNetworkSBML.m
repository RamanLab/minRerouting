function [  ] = SLNetworkSBML( model, Jdl, minRerouting, Clusters )
%SLNetworkSBML creates SBML file for reactions invovled in minimal
%rerouting with clusters as subsystems

% Unique of all minimal reroutings
SLNetRxn=[];
for iLeth=1:length(minRerouting)
    SLNetRxn=[SLNetRxn;minRerouting(iLeth).rxns];
end
SLNetRxn=unique(SLNetRxn);
SLNetRxnId=find(ismember(model.rxns,SLNetRxn));
SLNetMet=[];
for iRxn=1:length(SLNetRxnId)
    SLNetMet=[SLNetMet;model.mets(find(model.S(:,SLNetRxnId(iRxn))))];
end
SLNetMet=unique(SLNetMet);
SLNetMetId=find(ismember(model.mets,SLNetMet));

% Defining Subsystems
rxnClustMat=zeros(length(SLNetRxn),length(Clusters));

for iCluster=1:length(Clusters)
    minReSet=getMinSetofCluster(Clusters(iCluster).Jdl,Jdl,minRerouting);
    rxnClustMat(find(ismember(SLNetRxn,minReSet)),iCluster)=1;
end
S=sum(rxnClustMat,2);
% Reactions outside any cluster
rxnIn0ClustID=find(S==0);
rxnIn1ClustID=find(S==1);
rxnIn2ClustID=find(S==2);
rxnIn3ClustID=find(S==3);
% Genration of PseudoModel
model_SLNetwork=removeRxns(model,model.rxns(find(~ismember(model.rxns,SLNetRxn))));
model_SLNetwork.description=horzcat(model.description,'_SLNetwork');

% Subsystems replaced by clusterIDs
model_SLNetwork.subSystems=cell(length(SLNetRxnId),1);
model_SLNetwork.subSystems(rxnIn0ClustID)={'No_Cluster'};
for i=1:length(rxnIn1ClustID)
    model_SLNetwork.subSystems(rxnIn1ClustID(i))={sprintf('Cluster_%d',find(rxnClustMat(rxnIn1ClustID(i),:)))};
end
for i=1:length(rxnIn2ClustID)
    clustIDs=find(rxnClustMat(rxnIn2ClustID(i),:));
    % Change ID of first reaction assigning cluster
    model_SLNetwork.subSystems(rxnIn2ClustID(i))={sprintf('Cluster_%d',clustIDs(1))};
    model_SLNetwork.rxns(rxnIn2ClustID(i))={horzcat(model_SLNetwork.rxns{rxnIn2ClustID(i)},sprintf('_Cluster_%d',clustIDs(1)))};
    
    % Add extra reactions
    rxnName=horzcat(model_SLNetwork.rxns{rxnIn2ClustID(i)},sprintf('_Cluster_%d',clustIDs(2)));
    rxnFormula=printRxnFormula(model_SLNetwork,model_SLNetwork.rxns{rxnIn2ClustID(i)});
    subSystem={sprintf('Cluster_%d',clustIDs(2))};
    model_SLNetwork=addReaction(model_SLNetwork,rxnName,rxnFormula,'subSystem',subSystem)
end
end






