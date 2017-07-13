function [minRerouting] = minReroutingRxns(model,Jdl,cutOff)
    nLethals = length(Jdl);
    
    minRerouting(nLethals).rxns=[];
    minRerouting(nLethals).diff=[];
%     minRerouting(nLethals).rxns1=[];
%     minRerouting(nLethals).rxns2=[];
    minRerouting(nLethals).PathShort=[];
    minRerouting(nLethals).PathLong=[];
    minRerouting(nLethals).pathCommon=[];   
    
    for iLeth=1:nLethals
            rxn1=Jdl(iLeth,1); rxn2=Jdl(iLeth,2);
            [flux1,flux2] = myMOMA(model,rxn1,rxn2);
            diff=flux1-flux2;
            minRerouting(iLeth).rxns=model.rxns(find(abs(diff)>cutOff));
            minRerouting(iLeth).diff=diff(find(abs(diff)>cutOff));
            
            flux1Rxn=model.rxns(find(flux1));
            flux2Rxn=model.rxns(find(flux2));
            Path2=minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns,flux1Rxn));
            Path1=minRerouting(iLeth).rxns(ismember(minRerouting(iLeth).rxns,flux2Rxn));

            minRerouting(iLeth).pathCommon=Path1(ismember(Path1,Path2));
            path1_Ex=Path1(~ismember(Path1,flux1Rxn));
            path2_Ex=Path2(~ismember(Path2,flux2Rxn));
            
            if length(path1_Ex)<=length(path2_Ex)
                minRerouting(iLeth).PathShort=path1_Ex;
                minRerouting(iLeth).PathLong=path2_Ex;
            else
                minRerouting(iLeth).PathShort=path2_Ex;
                minRerouting(iLeth).PathLong=path1_Ex;
            end    
   
%             [flux1,flux2,fluxWT] = myMOMAvsWT(model,rxn1,rxn2);
%             diff1=flux1-fluxWT;
%             diff2=flux2-fluxWT;
%             
%             minRerouting(iLeth).rxns1=model.rxns(find(abs(diff1)>cutOff));
%             minRerouting(iLeth).rxns2=model.rxns(find(abs(diff2)>cutOff));           
           
    end
end
