### Subnetworks representing synthetic lethal submodules

1. For Ecoli iJO1366 model, 267 double lethal pairs were identified. 
2. 237 pairs among these were classified into 31 clusters based on common reactions in minimal rerouting using hierarchical clustering. 
3. Given common minimal rerouting set for a particular cluster, we identify minimal rerouting between all possible pairs in the set.
4. Classify the reaction pairs with no rerouting into one group and those with rerouting in different groups.

### Example - One of the clusters of 9 pairs have 17 reactions in rerouting 

|Lethal Pairs| |             
|--- | --- |                 
AACPS4	|3HAD161
AACPS4	|3OAR161
AACPS4	|3OAS161
CTECOAI7	|3HAD161
CTECOAI7	|3OAR161
CTECOAI7	|3OAS161
FACOAE161	|3HAD161
FACOAE161	|3OAR161
FACOAE161	|3OAS161

|Minimal Rerouting Set| 
|---|
3HAD161
3OAR161
3OAS161
AACPS2
AACPS4
ACACT7r
ACOAD6f
ACOATA
CTECOAI6
CTECOAI7
EAR161x
ECOAH7
FACOAE141
FACOAE161
FADRx2
HACD7
MACPD

We identify rerouting between each of the 136 pairs derived from set of 17 reactions.

| 3HAD161   | 3OAR161 | 3OAS161 | AACPS2 | AACPS4 | ACACT7r | ACOAD6f | ACOATA | CTECOAI6 | CTECOAI7 | EAR161x | ECOAH7 | FACOAE141 | FACOAE161 | FADRx2 | HACD7 | MACPD |   |
|------------|----------|----------|---------|---------|----------|----------|---------|-----------|-----------|----------|---------|------------|------------|---------|--------|--------|---|
| 3HAD161   | 0        | 0        | 0       | 0       | 17       | 17       | 17      | 0         | 0         | 17       | 0       | 17         | 0          | 17      | 0      | 17     | 0 |
| 3OAR161   | 0        | 0        | 0       | 0       | 17       | 17       | 17      | 0         | 0         | 17       | 0       | 17         | 0          | 17      | 0      | 17     | 0 |
| 3OAS161   | 0        | 0        | 0       | 0       | 17       | 17       | 17      | 0         | 0         | 17       | 0       | 17         | 0          | 17      | 0      | 17     | 0 |
| AACPS2    | 0        | 0        | 0       | 0       | 17       | 17       | 17      | 0         | 0         | 17       | 0       | 17         | 0          | 17      | 0      | 17     | 0 |
| AACPS4    | 17       | 17       | 17      | 17      | 0        | 0        | 0       | 0         | 17        | 0        | 0       | 0          | 17         | 0       | 0      | 0      | 0 |
| ACACT7r   | 17       | 17       | 17      | 17      | 0        | 0        | 0       | 0         | 17        | 0        | 0       | 0          | 17         | 0       | 0      | 0      | 0 |
| ACOAD6f   | 17       | 17       | 17      | 17      | 0        | 0        | 0       | 0         | 17        | 0        | 0       | 0          | 17         | 0       | 0      | 0      | 0 |
| ACOATA    | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0         | 0         | 0        | 0       | 0          | 0          | 0       | 0      | 0      | 0 |
| CTECOAI6  | 0        | 0        | 0       | 0       | 17       | 17       | 17      | 0         | 0         | 17       | 0       | 17         | 0          | 17      | 0      | 17     | 0 |
| CTECOAI7  | 17       | 17       | 17      | 17      | 0        | 0        | 0       | 0         | 17        | 0        | 0       | 0          | 17         | 0       | 0      | 0      | 0 |
| EAR161x   | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0         | 0         | 0        | 0       | 0          | 0          | 0       | 0      | 0      | 0 |
| ECOAH7    | 17       | 17       | 17      | 17      | 0        | 0        | 0       | 0         | 17        | 0        | 0       | 0          | 17         | 0       | 0      | 0      | 0 |
| FACOAE141 | 0        | 0        | 0       | 0       | 17       | 17       | 17      | 0         | 0         | 17       | 0       | 17         | 0          | 17      | 0      | 17     | 0 |
| FACOAE161 | 17       | 17       | 17      | 17      | 0        | 0        | 0       | 0         | 17        | 0        | 0       | 0          | 17         | 0       | 0      | 0      | 0 |
| FADRx2    | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0         | 0         | 0        | 0       | 0          | 0          | 0       | 0      | 0      | 0 |
| HACD7     | 17       | 17       | 17      | 17      | 0        | 0        | 0       | 0         | 17        | 0        | 0       | 0          | 17         | 0       | 0      | 0      | 0 |
| MACPD     | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0         | 0         | 0        | 0       | 0          | 0          | 0       | 0      | 0      | 0 |


Classifying subnetwork based on pairwise rerouting:

Reactions with zero rerouting with every other reaction in the subnetwork are those which belong to both the pathways.  
| Path 1     | Path 2     | Common   |
|------------|------------|----------|
| 3HAD161   | AACPS4    | ACOATA  |
| 3OAR161   | ACACT7r   | EAR161x |
| 3OAS161   | ACOAD6f   | FADRx2  |
| AACPS2    | CTECOAI7  | MACPD   |
| CTECOAI6  | ECOAH7    |         |
| FACOAE141 | FACOAE161 |         |
|           | HACD7     |         |
