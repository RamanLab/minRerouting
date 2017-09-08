### Subnetwork no 2 representing cluster no 2

| Lethal Pairs |         |
|--------------|---------|
| ALATA_L      | VPAMTr  |
| VALTA        | VPAMTr  |
| VALTA        | ALATA_L |


|Minimal Rerouting Set| 
|---|
'ALATA_L'
'VALTA'
'VPAMTr'

#### Rerouting between all the pairs
| Rxn      | ALATA_L' | VALTA' | VPAMTr' |
|----------|----------|--------|---------|
| ALATA_L' | 0        | 3      | 3       |
| VALTA'   | 3        | 0      | 3       |
| VPAMTr'  | 3        | 3      | 0       |

#### Classifying subnetwork based on pairwise rerouting:
| Path 1   | Path 2 | Path3   |
|----------|--------|---------|
| ALATA_L' | VALTA' | VPAMTr' |

- All are in parallel with each other
