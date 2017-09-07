### Subnetworks representing synthetic lethal submodules

1. For Ecoli iJO1366 model, 267 double lethal pairs were identified. 
2. 237 pairs among these were classified into 31 clusters based on common reactions in minimal rerouting using hierarchical clustering. 
3. Given common minimal rerouting set for a particular cluster, we identify minimal rerouting between all possible pairs in the set.
4. Classify the reaction pairs with no rerouting into one group and those with rerouting in different groups.

### Example - One of the clusters of 9 pairs have 17 reactions in rerouting 

|Lethal Pairs| |
|--- | ---|
AACPS4	|'3HAD161'
'AACPS4'	|'3OAR161'
'AACPS4'	|'3OAS161'
'CTECOAI7'	|'3HAD161'
'CTECOAI7'	|'3OAR161'
'CTECOAI7'	|'3OAS161'
'FACOAE161'	|'3HAD161'
'FACOAE161'	|'3OAR161'
'FACOAE161'	|'3OAS161'

We identify rerouting between each of the 136 pairs derived from set of 17 reactions.
