# minRerouting
Package to understand the compensatory mechanisms performing similar metabolic functions through analysis of Synthetic Lethals

# Alternate Routes in Networks
The rerouting of the reaction fluxes upon deletion of genes/reactions in the metabolic network is non-trivial and involves complex structrues. This package helps to identify minimal rerouting corresponding to synthetic lethal pairs using constraint based tools.   

[Workflow for analysis of synthetic lethals](docs/workflow.jpg)

# Steps in analysis
1. Get FastSL toolbox from https://github.com/RamanLab/FastSL

2. Run FastSL on COBRA model/s object to identify synthetic lethal pairs for given models.
Here we explain the example of some of the models of pathogenic organims downloaded from BIGG database. Exchange reactions are exluded from the analysis to focus on intracellular properties of metabolic networks

```Matlab 
>> addpath('core')
>> model_names = {'iIT341', 'iML1515', 'iNJ661', 'iPC815', 'iYL1228', 'STM_v1_0'} % Recommneded to rename model.mat file to modelNames.mat file
>> path_to_models = {'examples/iIT341/', 'examples/iML1515/', 'examples/iNJ661/','examples/iPC815/', 'examples/iYL1228/', 'examples/STM_v1_0/'}
>> Castle = getFastSL(model_names, path_to_models, 2); % 2 for double lethals and 3 for triple lethals
``` 

Castle.data has models, SLs, modelNames
 To Do: modify for without fastSL run in cases where SL are identified 

3. Run minRerouting
```Matlab 
>> Castle = multiMinRerouting(Castle, 'one')
``` 

4. Comparitive analysis synthetic lethals 
Use Castle input Struct to run compareLethals
```Matlab 
>> Castle = compareLethals(Castle, 2); 
```
Castle.data with new fields:
lethStatus 
lethCount 
all_rxns
core(Single, Double, Nonessential, Single_to_double)
simTable (SingleCnt, Single.rxns) 
 To Do: modify for TL compare

```Matlab 
>> [Castle, lethPair] = compareModelsPairwise(Castle, model1, model2)
```
To Do : Need to finish

5. Save the output data in xls
Use Castle structure to save output data in xls format
```Matlab 
>> get_output_castle(Castle_l1, 'examples/L1_minimal')
>> get_output_castle(Castle_l2, 'examples/L2_minimal')
```

