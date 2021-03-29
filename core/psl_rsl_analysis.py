import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

def combine_results(df1, df2):
    df2["Rxn_1_Type"] = df1["Rxn_1_Class"]
    df2["Rxn_2_Type"] = df1["Rxn_2_Class"]
    df2 = df2[["Rxn_1", "Rxn_1_Type", "Rxn_1_Min", "Rxn_1_Max", "Rxn_2", "Rxn_2_Type", "Rxn_2_Min", "Rxn_2_Max", "v1", "v2"]]
    pd.set_option("display.max_rows", None, "display.max_columns", None)
    
    return df2

def psl_rsl(df):
    df.loc[(np.abs(df["v1"] - 0)>1e-8) & (np.abs(df["v2"]-0)<1e-8), "Type"] = "PSL"
    df.loc[(np.abs(df["v1"] - 0)<1e-8) & (np.abs(df["v2"]-0)>1e-8), "Type"] = "PSL"
    df.loc[(np.abs(df["v1"] - 0)>1e-8) & (np.abs(df["v2"]-0)>1e-8), "Type"] = "RSL"
    
    return df

if __name__ == "__main__":
    df_list = []
    verbosity = 0

    for model in ['iIT341', 'iML1515', 'iNJ661', 'iPC815', 'iYL1228', 'STM_v1_0', 'e_coli_core']:
        print("PSL/RSL analysis for Organism", model, "... ", end="")
        df1 = pd.read_csv("../examples/" + model + "/" + model + "_pFBA.csv")
        df2 = pd.read_csv("../examples/" + model + "/" + model + "_FVA_one_norm.csv")
        df2 = combine_results(df1, df2)
        
        print("Done!")

        df2 = psl_rsl(df2)
        df_list.append(df2)
        df2.to_csv("../results/" + model + "/" + model + "_PSL_RSL_one_norm.csv")

        if verbosity: print(df2)