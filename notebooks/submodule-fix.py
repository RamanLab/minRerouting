import os
import numpy as np
import pandas as pd
from analysis.utils import get_data

# %load_ext autoreload
# %autoreload 2

df_new2_list = []
models = ['e_coli_core', 'iECW_1372', 'iIT341', 'iJO1366', 'iML1515', 'iNJ661', 'iPC815', 'iSSON_1240', 'iYL1228', 'iYS1720', 'STM_v1_0']

os.chdir("../results")
print(os.listdir())

for i in models:
    print("Organism:", i)
    df_new = pd.read_csv(i + "/" + i +'_submodules.csv')
    df_new2 = df_new.copy()
    df_new2.fillna("None", inplace=True)
    
    ss1 = df_new2["Subsystem 1"].to_list()
    ss2 = df_new2["Subsystem 2"].to_list()

    rxn1 = df_new2["Rxn_1"].to_list()
    rxn2 = df_new2["Rxn_2"].to_list()

    rxn1_new = []
    rxn2_new = []
    same_flag = []
    new = pd.DataFrame()

    pair_submodules = [(i, j) for i, j in zip(ss1, ss2)]
    for k, j in enumerate(pair_submodules):
        # print(j, tuple(sorted(j)))
        if j == tuple(sorted(j)):
            rxn1_new.append(rxn1[k])
            rxn2_new.append(rxn2[k])
            pair_submodules[k] = str(j[0] + " & " + j[1])
        else:
            rxn1_new.append(rxn2[k])
            rxn2_new.append(rxn1[k])
            pair_submodules[k] = str(j[1] + " & " + j[0])

        if j[0] == j[1]:
            same_flag.append(1)
        else:
            same_flag.append(0)

    new["Reaction 1"] = rxn1_new
    new["Reaction 2"] = rxn2_new
    new["Combination"] = pair_submodules
    new["Same"] = same_flag

    new = new.groupby(["Combination"]).count()
    new = new.reset_index()

    new["Count"] = new["Reaction 1"]
    new["Same"] = new["Combination"].apply(lambda x: x.split(" & ")[0]==x.split(" & ")[1])
    new = new[["Combination", "Count", "Same"]]
    new["Percentage"] = new["Count"]/np.sum(new["Count"])
    new.sort_values(["Count", "Percentage"], ascending=False, inplace=True)
    new.reset_index(drop=True, inplace=True)
    new.to_csv(i + "/" + i +"_analyzed_submodules.csv")

    value = new["Count"].sum()
    new = new.head()
    new_value = new["Count"].sum()
    print(value, new_value)

    temp_new = {"Combination":"All Other Submodule Combinations", "Count":value-new_value, "Percentage":1-(new_value/value)}
    new = new.append(temp_new, ignore_index=True)
    print(new)
    new.to_csv("csv/"+i+"_submodule.csv")
