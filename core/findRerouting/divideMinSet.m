function [Divide] = divideMinSet( model, minRerouting, Jdl )
%divideMinSet divides the minimal rerouting network in to alternate paths
%   Detailed explanation goes here
     Divide(length(Jdl)).Common=[];
     Divide(length(Jdl)).PathA_Exc=[];
     Divide(length(Jdl)).PathB_Exc=[];

    for iLeth=1:length(Jdl)
        if length(minRerouting(iLeth).rxns)~=0
            rxn1=find(strcmp(Jdl(iLeth,1),model.rxns));
            rxn2=find(strcmp(Jdl(iLeth,2),model.rxns));

            modeldel_1=model;
            modeldel_1.lb(rxn1)=0;
            modeldel_1.ub(rxn1)=0;

            modeldel_2=model;
            modeldel_2.lb(rxn2)=0;
            modeldel_2.ub(rxn2)=0;

            sol_1=optimizeCbModel(modeldel_1,'max','zero');
            sol_2=optimizeCbModel(modeldel_2,'max','zero'); 
            sol1rxn=model.rxns(find(sol_1.x));
            sol2rxn=model.rxns(find(sol_2.x));

            PathA=minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns,sol1rxn));
            PathB=minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns,sol2rxn));

            Divide(iLeth).Common=PathA(ismember(PathA,PathB));
            Divide(iLeth).PathA_Exc=PathA(~ismember(PathA,sol2rxn));
            Divide(iLeth).PathB_Exc=PathB(~ismember(PathB,sol1rxn));
        end 
    end   

end


