get_escher



% List of reactions in minRe universe (Sets less than 50)
uni_minRe = []
for iLeth = 1:length(minReSets)
    if length(minReSets(iLeth).rxns)<50
        uni_minRe = unique([uni_minRe;minReSets(iLeth).rxns]);
    end    
end
for k = 1:length(uni_minRe)
    if ismember(uni_minRe(k),Jdl(:,:))
        rank(k,1) = -1;
    else 
        rank(k,1) = 0;
    end
end    