### Subnetwork no 9 representing cluster no 9

| Lethal Pairs |         |
|--------------|---------|
| AACPS7       | 3HAD120 |
| AACPS7       | 3OAR120 |
| AACPS7       | 3OAS120 |
| FACOAE120    | 3HAD120 |
| FACOAE120    | 3OAR120 |
| FACOAE120    | 3OAS120 |

|Minimal Rerouting Set| 
|---|
'3HAD120'
'3OAR120'
'3OAS120'
'AACPS7'
'AACPS8'
'ACACT5r'
'ACOAD5f'
'ACOATA'
'EAR120x'
'ECOAH5'
'FACOAE100'
'FACOAE120'
'FADRx2'
'HACD5'
'MACPD'

#### Rerouting between all the pairs

| Rxns       | 3HAD120' | 3OAR120' | 3OAS120' | AACPS7' | AACPS8' | ACACT5r' | ACOAD5f' | ACOATA' | EAR120x' | ECOAH5' | FACOAE100' | FACOAE120' | FADRx2' | HACD5' | MACPD' |
|------------|----------|----------|----------|---------|---------|----------|----------|---------|----------|---------|------------|------------|---------|--------|--------|
| 3HAD120'   | 0        | 0        | 0        | 15      | 0       | 15       | 15       | 0       | 0        | 15      | 0          | 15         | 0       | 15     | 0      |
| 3OAR120'   | 0        | 0        | 0        | 15      | 0       | 15       | 15       | 0       | 0        | 15      | 0          | 15         | 0       | 15     | 0      |
| 3OAS120'   | 0        | 0        | 0        | 15      | 0       | 15       | 15       | 0       | 0        | 15      | 0          | 15         | 0       | 15     | 0      |
| AACPS7'    | 15       | 15       | 15       | 0       | 22      | 0        | 0        | 0       | 0        | 0       | 22         | 0          | 0       | 0      | 0      |
| AACPS8'    | 0        | 0        | 0        | 22      | 0       | 179      | 179      | 0       | 0        | 179     | 0          | 86         | 0       | 179    | 0      |
| ACACT5r'   | 15       | 15       | 15       | 0       | 179     | 0        | 0        | 0       | 0        | 0       | 22         | 0          | 0       | 0      | 0      |
| ACOAD5f'   | 15       | 15       | 15       | 0       | 179     | 0        | 0        | 0       | 0        | 0       | 22         | 0          | 0       | 0      | 0      |
| ACOATA'    | 0        | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0        | 0       | 0          | 0          | 0       | 0      | 0      |
| EAR120x'   | 0        | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0        | 0       | 0          | 0          | 0       | 0      | 0      |
| ECOAH5'    | 15       | 15       | 15       | 0       | 179     | 0        | 0        | 0       | 0        | 0       | 22         | 0          | 0       | 0      | 0      |
| FACOAE100' | 0        | 0        | 0        | 22      | 0       | 22       | 22       | 0       | 0        | 22      | 0          | 86         | 0       | 179    | 0      |
| FACOAE120' | 15       | 15       | 15       | 0       | 86      | 0        | 0        | 0       | 0        | 0       | 86         | 0          | 0       | 0      | 0      |
| FADRx2'    | 0        | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0        | 0       | 0          | 0          | 0       | 0      | 0      |
| HACD5'     | 15       | 15       | 15       | 0       | 179     | 0        | 0        | 0       | 0        | 0       | 179        | 0          | 0       | 0      | 0      |
| MACPD'     | 0        | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0        | 0       | 0          | 0          | 0       | 0      | 0      |


#### Classifying subnetwork based on pairwise rerouting:

| Path 1     | Path 2     | Common   | High Rewiring |
|------------|------------|----------|---------------|
| 3HAD120'   | AACPS7'    | ACOATA'  | AACPS8'       |
| 3OAR120'   | ACACT5r'   | EAR120x' |               |
| 3OAS120'   | ACOAD5f'   | FADRx2'  |               |
| FACOAE100' | ECOAH5'    | MACPD'   |               |
|            | FACOAE120' |          |               |
|            | HACD5'     |          |               |
