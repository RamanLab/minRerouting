function [ inputNetwork ] = getNetkwork( rxnSet,model )
%getNetwork returns the input file for Cytoscape visualization
% Input
% rxnSet - The set of reactions to be plotted in cyctoscape
% model - SBML model  
% Output
% inputNetwork - a matrix to be imported in csv format for Cytoscape
for k=1:length(rxnSet)
    rxnSetId(k)=find(strcmp(rxnSet(k),model.rxns));  
end    
rxnSetS=model.S(:,rxnSetId);
inputNetwork=[];
for m=1:length(rxnSet)
    if model.rev(rxnSetId(m))==1
        for p=1:length(model.mets)
            if rxnSetS(p,m)~=0
                inputNetwork=[inputNetwork;[model.mets(p),rxnSet(m),'species','reaction']];
                inputNetwork=[inputNetwork;[rxnSet(m),model.mets(p),'reaction','species']];
            end    
        end
    else
        for p=1:length(model.mets) 
            if rxnSetS(p,m)<0 
                inputNetwork=[inputNetwork;[model.mets(p),rxnSet(m),'species','reaction']];
            elseif rxnSetS(p,m)>0
                inputNetwork=[inputNetwork;[rxnSet(m),model.mets(p),'reaction','species']];
            end
        end  
    end 
end 

end

