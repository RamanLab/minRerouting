### Subnetwork no 11 representing cluster no 11

| Lethal Pairs |        |
|--------------|--------|
| PAPSR        | GRXR   |
| PAPSR        | GTHOr  |
| PAPSR        | PAPSR2 |
| TRDR         | GRXR   |
| TRDR         | GTHOr  |
| TRDR         | PAPSR2 |

|Minimal Rerouting Set| 
|---|
'GRXR'
'GTHOr'
'PAPSR'
'PAPSR2'
'TRDR'

#### Rerouting between all the pairs
| Rxn     | GRXR' | GTHOr' | PAPSR' | PAPSR2' | TRDR' |
|---------|-------|--------|--------|---------|-------|
| GRXR'   | 0     | 0      | 5      | 0       | 5     |
| GTHOr'  | 0     | 0      | 5      | 0       | 5     |
| PAPSR'  | 5     | 5      | 0      | 5       | 0     |
| PAPSR2' | 0     | 0      | 5      | 0       | 5     |
| TRDR'   | 5     | 5      | 0      | 5       | 0     |

#### Classifying subnetwork based on pairwise rerouting:
| Path 1  | Path 2 | Common |
|---------|--------|--------|
| GRXR'   | PAPSR' |        |
| GTHOr'  | TRDR'  |        |
| PAPSR2' |        |        |
