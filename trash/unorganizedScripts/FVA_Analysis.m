initCobraToolbox; 
load('iAF1260Cplex_Rxn_lethals');
load('iAF1260');
model=iAF1260;

Jlr=unique(Jdl(:,:));
nLethals = length(Jlr);

% FVA on WT
[minFluxWT, maxFluxWT]= fluxVariability(model,5);

% FVA on unique reactions in lethal pair
fva_Leth_5(nLethals).minFlux=[];
fva_Leth_5(nLethals).maxFlux=[];
h = waitbar(0,'0.00','Name','FVA on lethals...'); 

for iLeth=1:nLethals
    rxn=find(strcmp(Jlr(iLeth),model.rxns));
   
    modeldel=model;
    modeldel.lb(rxn)=0;
    modeldel.ub(rxn)=0;
    
    [fva_Leth_5(iLeth).minFlux,fva_Leth_5(iLeth).maxFlux] = fluxVariability(modeldel,5);
    waitbar(iLeth/length(Jlr),h,[num2str(round(iLeth*100/length(Jlr))) '% completed...']);

end

minRerouting_FVA(length(Jdl)).rxns=[];
for iLeth=1:length(Jdl)
    rxn1=find(strcmp(Jdl(iLeth,1),Jlr));
    rxn2=find(strcmp(Jdl(iLeth,2),Jlr));
    % FVA on MT1 and MT2
    [minFluxMT1, maxFluxMT1]= deal(fva_Leth_5(rxn1).minFlux,fva_Leth_5(rxn1).maxFlux);
    [minFluxMT2, maxFluxMT2]= deal(fva_Leth_5(rxn2).minFlux,fva_Leth_5(rxn2).maxFlux);
    % Overlap
    LB=max(minFluxMT1,minFluxMT2);
    UB=min(maxFluxMT1,maxFluxMT2);
    overlap=UB-LB;
    range=max(maxFluxMT1,maxFluxMT2)-min(minFluxMT1,minFluxMT2);

    OF=overlap./range;

%     % Rxns with range=0
%     zeroRange=model.rxns(range==0);
%     zeroOverlap=model.rxns(overlap<=0);
%     % Rxns with nonzero range but zero overlap
%     minSet=zeroOverlap(find(~ismember(zeroOverlap,zeroRange)));
    minRerouting_FVA(iLeth).rxns=model.rxns(find(OF==0))
end    