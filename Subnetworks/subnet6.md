### Subnetwork no 2 representing cluster no 2
| Lethal Pairs |         |
|--------------|---------|
| ACACT4r      | 3HAD100 |
| ACACT4r      | 3OAR100 |
| ACACT4r      | 3OAS100 |
| ACOAD4f      | 3HAD100 |
| ACOAD4f      | 3OAR100 |
| ACOAD4f      | 3OAS100 |
| ECOAH4       | 3HAD100 |
| ECOAH4       | 3OAR100 |
| ECOAH4       | 3OAS100 |
| HACD4        | 3HAD100 |
| HACD4        | 3OAR100 |
| HACD4        | 3OAS100 |

|Minimal Rerouting Set| 
|---|
'3HAD100'
'3OAR100'
'3OAS100'
'AACPS8'
'AACPS9'
'ACACT4r'
'ACOAD4f'
'ACOATA'
'EAR100x'
'ECOAH4'
'FACOAE100'
'FACOAE80'
'FADRx2'
'HACD4'
'MACPD'

#### Rerouting between all the pairs
| | 3HAD100'   | 3OAR100' | 3OAS100' | AACPS8' | AACPS9' | ACACT4r' | ACOAD4f' | ACOATA' | EAR100x' | ECOAH4' | FACOAE100' | FACOAE80' | FADRx2' | HACD4' | MACPD' |  
|------------|----------|----------|---------|---------|----------|----------|---------|----------|---------|------------|-----------|---------|--------|--------|---|
| 3HAD100'   | 0        | 0        | 0       | 0       | 0        | 15       | 15      | 0        | 0       | 15         | 0         | 0       | 0      | 15     | 0 |
| 3OAR100'   | 0        | 0        | 0       | 0       | 0        | 15       | 15      | 0        | 0       | 15         | 0         | 0       | 0      | 15     | 0 |
| 3OAS100'   | 0        | 0        | 0       | 0       | 0        | 15       | 15      | 0        | 0       | 15         | 0         | 0       | 0      | 15     | 0 |
| AACPS8'    | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0        | 0       | 0          | 0         | 0       | 0      | 0      | 0 |
| AACPS9'    | 0        | 0        | 0       | 0       | 0        | 15       | 15      | 0        | 0       | 15         | 0         | 0       | 0      | 15     | 0 |
| ACACT4r'   | 15       | 15       | 15      | 0       | 15       | 0        | 0       | 0        | 0       | 0          | 0         | 15      | 0      | 0      | 0 |
| ACOAD4f'   | 15       | 15       | 15      | 0       | 15       | 0        | 0       | 0        | 0       | 0          | 0         | 15      | 0      | 0      | 0 |
| ACOATA'    | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0        | 0       | 0          | 0         | 0       | 0      | 0      | 0 |
| EAR100x'   | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0        | 0       | 0          | 0         | 0       | 0      | 0      | 0 |
| ECOAH4'    | 15       | 15       | 15      | 0       | 15       | 0        | 0       | 0        | 0       | 0          | 0         | 15      | 0      | 0      | 0 |
| FACOAE100' | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0        | 0       | 0          | 0         | 0       | 0      | 0      | 0 |
| FACOAE80'  | 0        | 0        | 0       | 0       | 0        | 15       | 15      | 0        | 0       | 15         | 0         | 0       | 0      | 15     | 0 |
| FADRx2'    | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0        | 0       | 0          | 0         | 0       | 0      | 0      | 0 |
| HACD4'     | 15       | 15       | 15      | 0       | 15       | 0        | 0       | 0        | 0       | 0          | 0         | 15      | 0      | 0      | 0 |
| MACPD'     | 0        | 0        | 0       | 0       | 0        | 0        | 0       | 0        | 0       | 0          | 0         | 0       | 0      | 0      | 0 |

#### Classifying subnetwork based on pairwise rerouting:

| Path 1   | Path 2    | Common     |
|----------|-----------|------------|
| 3HAD100' | ACACT4r'  | AACPS8'    |
| 3OAR100' | ACOAD4f'  | ACOATA'    |
| 3OAS100' | ECOAH4'   | EAR100x'   |
| AACPS9'  | FACOAE80' | FACOAE100' |
|          | HACD4'    | FADRx2'    |
|          |           | MACPD'     |
