import os
import numpy as np
import pandas as pd
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def combine_results(df_pFBA, df_FVA):
    """Function to combine the results obtained form pFBA and FVA."""

    # Obtain the reaction types
    df_FVA["Rxn_1_Type"] = df_pFBA["Rxn_1_Class"]
    df_FVA["Rxn_2_Type"] = df_pFBA["Rxn_2_Class"]

    # Rearrange the order of columns
    df_FVA = df_FVA[["Rxn_1", "Rxn_1_Type", "Rxn_1_Min", "Rxn_1_Max", "Rxn_2", "Rxn_2_Type", "Rxn_2_Min", "Rxn_2_Max", "v1", "v2"]]
    # Ensure that the whole dataframe can be printed
    pd.set_option("display.max_rows", None, "display.max_columns", None, "display.expand_frame_repr", False)
    
    return df_FVA

def psl_rsl(df):
    """Function to perform PSL/RSL analysis"""
    mask = (np.abs(df_FVA["Rxn_1_Min"])>1e-8) & (np.abs(df_FVA["Rxn_2_Min"])>1e-8)
    df.loc[mask, "Type"] = "RSL"
    df.loc[~(mask), "Type"] = "PSL"
    
    return df

def get_distribution(df, model):
    """Function to group the dataframe by the reaction types and the pair type"""
    df_grouped = df.groupby(["Type", "Rxn_1_Type", "Rxn_2_Type"]).count()
    df_grouped = df_grouped.rename(columns={"Rxn_1":"Count", "Rxn_1_Min":"Fraction"})
    df_grouped["Fraction"] /= df_grouped["Fraction"].sum()
    df_grouped_list.append(df_grouped)
    df_grouped.drop(["Rxn_1_Max", "Rxn_2", "Rxn_2_Min", "Rxn_2_Max", "v1", "v2"], axis=1, inplace=True)
    df_grouped.to_csv("../results/"+model+"/"+model+"_PSL_RSL_distribution.csv")
    df_grouped.reset_index(inplace=True)
    return df_grouped

def plot_results(df, model):
    """Function to plot the results"""
    df["Combination"] = df["Rxn_1_Type"] + ", " + df["Rxn_2_Type"]

    psl_df = df.loc[df["Type"]=="PSL",:]
    rsl_df = df.loc[df["Type"]=="RSL",:]

    fig = make_subplots(rows=1, cols=2, subplot_titles=("PSL", "RSL"), specs=[[{"type": "domain"}, {"type": "domain"}]])
    fig.add_trace(go.Pie(labels=psl_df["Combination"], values=psl_df["Fraction"]*100, hole=0.4),row=1, col=1)
    fig.add_trace(go.Pie(labels=rsl_df["Combination"], values=rsl_df["Fraction"]*100, hole=0.4),row=1, col=2)

    fig.update_layout(title_text="Model: "+model)
    image_name = "../results/images/"+model+"_PSL_RSL_dist.png"
    fig.write_image(image_name)
    fig.show()

if __name__ == "__main__":
    df_list = []
    verbosity = 0
    df_grouped_list = []

    for model in ['iIT341', 'iML1515', 'iNJ661', 'iPC815', 'iYL1228', 'STM_v1_0', 'e_coli_core']:
        print("PSL/RSL analysis for Organism", model, "... ", end="")
        # Obtain pFBA class of the reaction pair
        df_pFBA = pd.read_csv("../examples/" + model + "/" + model + "_pFBA.csv")
        # Get the flux distributions
        df_FVA = pd.read_csv("../examples/" + model + "/" + model + "_FVA_one_norm.csv")
        # Combine the results into one dataframe
        df_FVA = combine_results(df_pFBA, df_FVA)
        

        # Perform PSL/RSL segregation
        df_FVA = psl_rsl(df_FVA)
        df_list.append(df_FVA)
        df_FVA.to_csv("../results/" + model + "/" + model + "_PSL_RSL_one_norm.csv")

        # Get the distribution of the reaction pairs
        df_grouped = get_distribution(df_FVA, model)
        df_grouped_list.append(df_grouped)
        plot_results(df_grouped, model)
        print("Done!")

        if verbosity: 
            print("Distribution:")
            print(df_grouped)
            print("="*50)