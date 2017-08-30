## Comparing synthetic lethals classification by Güell et al (2014)

This excercise attempts to analyse and comparer the classification of lethals as proposed by [Güell et al (2014)](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1003637#s4). The hypothesis is that synthetic lethals can be classified into plasticity and redundancy subtypes in metabolic networks. 

### Classification:
* Plasticity SL (PSL) pairs, with one reaction functional and one silent in the WT
* Redundancy SL (RSL) pairs, with both reactions active

![Classification Lethals](https://github.com/RamanLab/minRerouting/blob/master/Comparison/Images/10.1371%252Fjournal.pcbi.1003637.g001.png)

This classification was based on the optimal FBA solution to maximize biomass growth rate. Metabolic constraints define a range of solution space. FBA is a based on linear programming approach to identify a particular solution maximizing given objective function, generally biomas  gorwth rate. Even though the solution space is reduced, FBA picks one of the solutions satisfying optimal growth. `optimizeCbModel` function in COBRA toolbox can be used to find such FBA solutions with additional objectives such as `minNorm`. We here compare several FBA solutions found using several minNorm conditions as well as different solvers.

#### Synthetic Lethals Pairs
We enlist synthetic lethal pairs identified using **Fast-SL** and as reported by Güell et al (2014).

| Organism      | Fast-SL       | Güell et al (2014) |
| :-------------: |:-----------:| :-----:|
| Ecoli iJO1366 | 267           | 256 |

List of additional lethal pairs identified using Fast-SL:

| Extra Pairs by Fast-SL  |
|:--- | :---|
| O2tex | OPHHX3  |
|O2tpp	|CPPPGO2|
|O2tpp	|OPHHX3|
|O2tpp	|SUCCtex|
|OPHHX	|OPHHX3|
|PDX5PO2	|PDX5POi|
|PRPPS	|R1PK|
|O2tpp	|ENO|
|O2tpp	|GAPD|
|PGK	|O2tpp|
|PGM	|O2tpp|
|PPPGO3	|O2tpp |


