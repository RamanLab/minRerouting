import os
import requests
import numpy as np
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import plotly.express as px
from scipy.io import loadmat
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from bs4 import BeautifulSoup as bs
from plotly.subplots import make_subplots

# Ensure that the whole dataframe can be printed
pd.set_option("display.max_rows", None, "display.max_columns", None, "display.expand_frame_repr", False)

import warnings
warnings.filterwarnings("ignore")
# plt.rcParams['font.size'] = 14
plt.rcParams['font.serif'] = "Cambria"
plt.rcParams['font.family'] = "serif"


def psl_rsl(df, model, tol=0):
    # Function to perform PSL/RSL analysis.
    # Ensure that both the minimum fluxes are greater than 0
    df["Rxn_1_Prod"] = df["Rxn_1_Min"]*df["Rxn_1_Max"]
    df["Rxn_2_Prod"] = df["Rxn_2_Min"]*df["Rxn_2_Max"]

    rsl_mask = (df["Rxn_1_Prod"] > 0) & (df["Rxn_2_Prod"] > 0)
    df.loc[rsl_mask, "Type"] = "RSL"

    df_ambi = pd.read_csv("../results/"+model+"/"+model+"_ambiguous_FVA_100.csv")
    df_ambi["g_prod"] = df_ambi["g_min_Fluxes"]*df_ambi["g_max_Fluxes"]
    df_ambi["l_prod"] = df_ambi["l_min_Fluxes"]*df_ambi["l_max_Fluxes"]

    psl_mask = (df_ambi["g_prod"] <= 0) | (df_ambi["l_prod"] <= 0)
    df_select = df.loc[~rsl_mask, :]
    df_select.reset_index(inplace=True)
    df_select.loc[psl_mask, "Type"] = "PSL"
    df_select.loc[~psl_mask, "Type"] = "RSL"
    df_select.set_index("index", inplace=True)

    df_new = pd.DataFrame()
    df_new = df.loc[rsl_mask]
    df_select = df_select.append(df_new)
    df_select.sort_index(inplace=True)

    return df_select

def combine_results(df_pFBA, df_FVA):
    '''
    Function to combine the results obtained 
    from pFBA and FVA.
    '''

    # Obtain the reaction types
    df_FVA["Rxn_1_Type"] = df_pFBA["Rxn_1_Class"]
    df_FVA["Rxn_2_Type"] = df_pFBA["Rxn_2_Class"]

    # Rearrange the order of columns
    df_FVA = df_FVA[["Rxn_1", "Rxn_1_Type", "Rxn_1_Min", \
                     "Rxn_1_Max", "Rxn_2", "Rxn_2_Type", \
                     "Rxn_2_Min", "Rxn_2_Max", "v1", "v2"]]
    
    return df_FVA

def get_distribution(df, model):
    # Function to group the dataframe by the reaction types and the pair type
    df_grouped = df.groupby(["Type", "Rxn_1_Type", "Rxn_2_Type"]).count()
    df_grouped = df_grouped.rename(columns={"Rxn_1":"Count", "Rxn_1_Min":"Fraction"})
    df_grouped["Fraction"] /= df_grouped["Fraction"].sum()
    df_grouped.drop(["Rxn_1_Max", "Rxn_2", "Rxn_2_Min", "Rxn_2_Max", "v1", "v2"], axis=1, inplace=True)
    df_grouped.to_csv("../results/"+model+"/"+model+"_PSL_RSL_distribution.csv")
    df_grouped.reset_index(inplace=True)
    return df_grouped

def plot_results(df_full, model_list):
    # Function to plot the results
    # Reference: https://stackoverflow.com/questions/8389636/creating-over-20-unique-legend-colors-using-matplotlib
    df_full["Combination"] = df_full["Rxn_1_Type"] + ", " + df_full["Rxn_2_Type"]

    unique_combination = list(df_full["Combination"].unique())
    unique_combination.sort()

    NUM_COLORS = len(unique_combination)
    sns.reset_orig()  
    clrs = sns.color_palette('husl', n_colors=NUM_COLORS)  # a list of RGB tuples
    clrs = np.array(clrs)
    clrs *= 255
    clrs = clrs.astype(int)
    clrs = list(clrs)
    clrs = [tuple(i) for i in clrs]
    clrs = ['rgb'+str(i) for i in clrs]

    color_mapping = {}
    for i,j in enumerate(unique_combination):
        color_mapping[j] = clrs[i]

    df_full["Colors"] = df_full["Combination"].map(color_mapping)

    count_list = []
    fraction_list = []

    for model in model_list:
        df = df_full.loc[df_full["Organism"]==model]

        psl_df = df.loc[df["Type"]=="PSL",:]
        rsl_df = df.loc[df["Type"]=="RSL",:]
        
        fig = make_subplots(rows=1, cols=2, subplot_titles=("PSL", "RSL"), \
                            specs=[[{"type": "domain"}, {"type": "domain"}]])
        
        fig.add_trace(go.Pie(labels=psl_df["Combination"], values=psl_df["Fraction"]*100, \
                             hole=0.4, marker_colors=psl_df["Colors"]), row=1, col=1)
        fig.add_trace(go.Pie(labels=rsl_df["Combination"], values=rsl_df["Fraction"]*100, \
                             hole=0.4, marker_colors=rsl_df["Colors"]), row=1, col=2)

        fig.update_layout(title_text="Model: "+model)
        image_name = "../results/images/"+model+"_PSL_RSL_dist.pdf"
        
        try:
            path = "../results/images/"
            try:
                os.mkdir(path)
            except:
                pass
            fin = open(image_name, "w")
            fin.close()
            fig.write_image(image_name)

        except:
            path = "../results/images/"
            try:
                os.mkdir(path)
            except:
                pass
            fig.write_image(image_name)
        
        fig.show()
        
        df_new = df.groupby(["Type"]).sum()
        # Ensure that order is PSL, RSL
        df_new.sort_index(inplace=True)
        count = df_new["Count"]
        fraction = df_new["Fraction"]
    
        count_list.append(count)
        fraction_list.append(fraction)

    return count_list, fraction_list

def get_diff(val):
    list_value = []
    for i in val:
        list_value.append(np.sum(np.sign(i)))
    
    return list_value

def get_net_diff(val):
    list_value = []
    for i in val:
        list_value.append(np.sum(i))
        
    return list_value


def preprocess(data):
    df = pd.DataFrame(data, columns=["rxns", "diff_flux", "abs_diff_flux", \
                                     "PathShort", "PathLong", "pathCommon", \
                                     "del_rxn1", "del_rxn2", "totalFluxDiff",\
                                     "solStatus"])

    # Get the reactions
    rxns = df["rxns"]
    new_rxns = []
    for i in rxns:
        temp = []
        for j in i:
            temp.append(str(j[0][0]))
        new_rxns.append(temp)
        
    df["rxns"] = new_rxns

    # Get reaction 1
    del_rxn1 = df["del_rxn1"]
    new_del_rxn1 = []
    for i in del_rxn1:
        new_del_rxn1.append(i[0][0][0])
    df["del_rxn1"] = new_del_rxn1

    # Get reaction 2
    del_rxn2 = df["del_rxn2"]
    new_del_rxn2 = []
    for i in del_rxn2:
        new_del_rxn2.append(i[0][0][0])
    df["del_rxn2"] = new_del_rxn2

    # Get flux difference
    diff_flux = df["diff_flux"]
    new_diff_flux = []
    for i in diff_flux:
        temp = []
        for j in i:
            temp.append(j[0])
        new_diff_flux.append(temp)
        
    df["diff_flux"] = new_diff_flux

    # Get absolute flux difference
    abs_diff_flux = df["abs_diff_flux"]
    new_abs_diff_flux = []
    for i in abs_diff_flux:
        temp = []
        for j in i:
            temp.append(j[0])
        new_abs_diff_flux.append(temp)
        
    df["abs_diff_flux"] = new_abs_diff_flux

    # Get PathShort
    PathShort = df["PathShort"]
    new_PathShort = []
    for i in PathShort:
        temp = []
        for j in i:
            temp.append(str(j[0][0]))
        new_PathShort.append(temp)
        
    df["PathShort"] = new_PathShort
    
    # Get PathLong
    PathLong = df["PathLong"]
    new_PathLong = []
    for i in PathLong:
        temp = []
        for j in i:
            temp.append(str(j[0][0]))
        new_PathLong.append(temp)
        
    df["PathLong"] = new_PathLong
    
    # Get pathCommon
    pathCommon = df["pathCommon"]
    new_pathCommon = []
    for i in pathCommon:
        temp = []
        for j in i:
            temp.append(str(j[0][0]))
        new_pathCommon.append(temp)
        
    df["pathCommon"] = new_pathCommon
    
    # Get totalFluxDiff
    totalFluxDiff = df["totalFluxDiff"]
    new_totalFluxDiff = []
    for i in totalFluxDiff:
        new_totalFluxDiff.append(i[0][0])
        
    df["totalFluxDiff"] = new_totalFluxDiff

    # Get solStatus
    solStatus = df["solStatus"]
    new_solStatus = []
    for i in solStatus:
        new_solStatus.append(i[0][0])
        
    df["solStatus"] = new_solStatus

    sl_size = [len(i) for i in df["rxns"]]
    df["sl_size"] = sl_size

    common_sl_size = [len(i) for i in df["pathCommon"]]
    df["common_sl_size"] = common_sl_size
    df["num_diff"] = get_diff(df["diff_flux"])
    df["net_diff"] = get_net_diff(df["diff_flux"])
    
    return df

def consolidate_results(model_list, norm):
    psl_df_main = pd.DataFrame()
    df_main = pd.DataFrame()
    # print(model_list)

    for model in tqdm(model_list):
        mat = loadmat("../results/" + model + "/Castle_" + norm + "_norm.mat")
        # Get the minRe data
        data = mat["data"][0][0][-1][0]
        df = preprocess(data)
        df["Organism"] = model
        
        df_new = pd.read_csv("../results/" + model + "/" + model + "_PSL_RSL_one_norm_100.csv", index_col=0)
        df["Type"] = df_new["Type"]
        df["Rxn_1_Active"] = np.abs(df_new["Rxn_1_Min"]) > 0
        df["Rxn_2_Active"] = np.abs(df_new["Rxn_2_Min"]) > 0
        df["Rxn_1_Min"] = df_new["Rxn_1_Min"]
        df["Rxn_2_Min"] = df_new["Rxn_2_Min"]
        df["v1"] = df_new["v1"]
        df["v2"] = df_new["v2"]
        
        # Same for psl_rsl dataframe
        psl_df = df.loc[df["Type"]=="PSL",:].copy()
        psl_rxn1_mask = (psl_df["Rxn_1_Active"] == True)   
        psl_df.loc[psl_rxn1_mask, "num_diff"] = psl_df.loc[psl_rxn1_mask, "num_diff"]*-1
        psl_df.loc[psl_rxn1_mask, "net_diff"] = psl_df.loc[psl_rxn1_mask, "net_diff"]*-1
        psl_df_main = psl_df_main.append(psl_df)
        
        # Get the PSL, rxn1 active and invert sign when rxn1 is active    
        # Because diff = flux1 - flux2 (< 0 => flux2 > flux 1)
        df_psl_rxn1 = ((df["Type"]=="PSL")&(df["Rxn_1_Active"]==True))
        df.loc[df_psl_rxn1, "num_diff"] = df.loc[df_psl_rxn1, "num_diff"]*-1
        df.loc[df_psl_rxn1, "net_diff"] = df.loc[df_psl_rxn1, "net_diff"]*-1
        
        # psl_df["Active_Reaction"] = np.select([rxn1_mask, rxn2_mask], [psl_df["Rxn_1"], psl_df["Rxn_2"]], default=np.nan)    
        # grouped_psl = psl_df.groupby(["Active_Reaction"])["Active_Reaction"].count()
        # psl_counts[model] = list(grouped_psl)
        
        df_select = df[["Organism", "sl_size", "common_sl_size", "num_diff", "net_diff"]]
        df_main = df_main.append(df_select)
        # display(df_new)
        # display(psl_df)

        fname = "../results/" + model + "/consolidated_" + norm + "_norm.csv"
        df.to_csv(fname)

        fname = "../results/" + model + "/analysis_" + norm + "_norm.csv"
        df_main.to_csv(fname)

        fname = "../results/" + model + "/psl_analysis_" + norm + "_norm.csv"
        psl_df_main.to_csv(fname)


    return df_main, df_select, psl_df_main

def fba_consolidate_results(model_list, norm):
    psl_df_main = pd.DataFrame()
    df_main = pd.DataFrame()
    # print(model_list)

    for model in tqdm(model_list):
        mat = loadmat("../results/" + model + "/fba_" + norm + "_norm.mat")
        # Get the minRe data
        data = mat["data"][0][0][-1][0]
        df = fba_preprocess(data)
        df["Organism"] = model
        
        df_new = pd.read_csv("../results/" + model + "/" + model + "_PSL_RSL_one_norm_100.csv", index_col=0)
        df["Type"] = df_new["Type"]
        df["Rxn_1_Active"] = np.abs(df_new["Rxn_1_Min"]) > 0
        df["Rxn_2_Active"] = np.abs(df_new["Rxn_2_Min"]) > 0
        df["Rxn_1_Min"] = df_new["Rxn_1_Min"]
        df["Rxn_2_Min"] = df_new["Rxn_2_Min"]
        df["v1"] = df_new["v1"]
        df["v2"] = df_new["v2"]
        
        # Same for psl_rsl dataframe
        psl_df = df.loc[df["Type"]=="PSL",:].copy()
        psl_rxn1_mask = (psl_df["Rxn_1_Active"] == True)   
        psl_df.loc[psl_rxn1_mask, "num_diff"] = psl_df.loc[psl_rxn1_mask, "num_diff"]*-1
        psl_df.loc[psl_rxn1_mask, "net_diff"] = psl_df.loc[psl_rxn1_mask, "net_diff"]*-1
        psl_df_main = psl_df_main.append(psl_df)
        
        # Get the PSL, rxn1 active and invert sign when rxn1 is active    
        # Because diff = flux1 - flux2 (< 0 => flux2 > flux 1)
        df_psl_rxn1 = ((df["Type"]=="PSL")&(df["Rxn_1_Active"]==True))
        df.loc[df_psl_rxn1, "num_diff"] = df.loc[df_psl_rxn1, "num_diff"]*-1
        df.loc[df_psl_rxn1, "net_diff"] = df.loc[df_psl_rxn1, "net_diff"]*-1
        
        # psl_df["Active_Reaction"] = np.select([rxn1_mask, rxn2_mask], [psl_df["Rxn_1"], psl_df["Rxn_2"]], default=np.nan)    
        # grouped_psl = psl_df.groupby(["Active_Reaction"])["Active_Reaction"].count()
        # psl_counts[model] = list(grouped_psl)
        
        df_select = df[["Organism", "sl_size", "common_sl_size", "num_diff", "net_diff"]]
        df_main = df_main.append(df_select)

        fname = "../results/" + model + "/fba_consolidated_" + norm + "_norm.csv"
        df.to_csv(fname)

        fname = "../results/" + model + "/fba_analysis_" + norm + "_norm.csv"
        df_main.to_csv(fname)


        fname = "../results/" + model + "/psl_fba_analysis_" + norm + "_norm.csv"
        psl_df_main.to_csv(fname)


    return df_main, df_select, psl_df_main

def consolidate_results_print_dfs(df_main, norm):
    from matplotlib.lines import Line2D
    line1 = Line2D([0], [0], label='Mean', color='k')
    line2 = Line2D([0], [0], label='Median', color='k', ls="--")

    columns_of_interest = ["sl_size", "common_sl_size", "num_diff", "net_diff"]
    fname_list = ["sl_cluster_size", "common_sl_size", "num_diff", "net_diff"]
    title_list = ["Distribution of SL cluster size", "Distribution of Common SL cluster size", \
                  "Distribution of number of fluxes with difference", "Distribution of total flux increase"]
    xlabel_list = ["# Reactions", "# Reactions", "# Reactions with different fluxes", "Total flux increase"]

    for col, fname, title, xlabel in zip(columns_of_interest, fname_list, title_list, xlabel_list):  
        df_select = df_main[[col, "Organism"]]

        plt.figure(figsize=(20,10))
        ax = sns.violinplot(x=col, y="Organism", data=df_main, color="0.8")
                            # palette=sns.color_palette("bright"), 
        ax = sns.stripplot(x=col, y="Organism", data=df_main, 
                        palette=sns.color_palette("bright"), jitter=True)

        ax = sns.pointplot(x=col, y='Organism', ci=None, data=df_main, 
                        color="k", estimator=np.mean, label="mean", scale=0.5) 
        ax = sns.pointplot(x=col, y='Organism', ci=None, data=df_main, color="k", 
                        linestyles='--', estimator=np.median, label="mean", scale=0.5) 

        plt.grid()
        ax.set_title(title)
        ax.set_ylabel("Organisms")
        ax.set_xlabel(xlabel)
        plt.legend(handles=[line1, line2])
        plt.savefig("../results/images/" + norm + "_norm_violin_"+fname+".png")
        plt.show()

        mean_df = df_main[[col, "Organism"]].groupby("Organism").mean()
        median_df = df_main[[col, "Organism"]].groupby("Organism").median()

        consolidated_df = mean_df.copy()
        mean_name = "mean_" + col
        median_name = "median_" + col
        consolidated_df[median_name] = median_df[col]
        consolidated_df.rename(columns={col:mean_name})

        try:
            consolidated_df.to_csv("../results/csv/" + norm + "_norm_" + fname + ".csv")
        except:
            os.mkdir("../results/csv/")
            consolidated_df.to_csv("../results/csv/" + norm + "_norm_" + fname + ".csv")
        display(consolidated_df)


def fba_consolidate_results_print_dfs(df_main, norm):
    from matplotlib.lines import Line2D
    line1 = Line2D([0], [0], label='Mean', color='k')
    line2 = Line2D([0], [0], label='Median', color='k', ls="--")

    columns_of_interest = ["sl_size", "common_sl_size", "num_diff", "net_diff"]
    fname_list = ["fba_sl_cluster_size", "fba_common_sl_size", "fba_num_diff", "fba_net_diff"]
    title_list = ["FBA Distribution of SL cluster size", \
                  "FBA Distribution of Common SL cluster size", \
                  "FBA Distribution of number of fluxes with difference", \
                  "FBA Distribution of total flux increase"]
    xlabel_list = ["# Reactions", "# Reactions", "# Reactions with different fluxes", "Total flux increase"]

    for col, fname, title, xlabel in zip(columns_of_interest, fname_list, title_list, xlabel_list):  
        df_select = df_main[[col, "Organism"]]

        plt.figure(figsize=(20,10))
        ax = sns.violinplot(x=col, y="Organism", data=df_main, color="0.8")
                            # palette=sns.color_palette("bright"), 
        ax = sns.stripplot(x=col, y="Organism", data=df_main, 
                        palette=sns.color_palette("bright"), jitter=True)

        ax = sns.pointplot(x=col, y='Organism', ci=None, data=df_main, 
                        color="k", estimator=np.mean, label="mean", scale=0.5) 
        ax = sns.pointplot(x=col, y='Organism', ci=None, data=df_main, color="k", 
                        linestyles='--', estimator=np.median, label="mean", scale=0.5) 

        plt.grid()
        ax.set_title(title)
        ax.set_ylabel("Organisms")
        ax.set_xlabel(xlabel)
        plt.legend(handles=[line1, line2])
        plt.savefig("../results/images/" + norm + "_norm_violin_"+fname+".png")
        plt.show()

        mean_df = df_main[[col, "Organism"]].groupby("Organism").mean()
        median_df = df_main[[col, "Organism"]].groupby("Organism").median()

        consolidated_df = mean_df.copy()
        mean_name = "mean_" + col
        median_name = "median_" + col
        consolidated_df[median_name] = median_df[col]
        consolidated_df.rename(columns={col:mean_name})

        try:
            consolidated_df.to_csv("../results/csv/" + norm + "_norm_" + fname + ".csv")
        except:
            os.mkdir("../results/csv/")
            consolidated_df.to_csv("../results/csv/" + norm + "_norm_" + fname + ".csv")
        display(consolidated_df)

def fba_preprocess(data):
    df = pd.DataFrame(data, columns=["rxns", "diff_flux", "abs_diff_flux", \
                                     "PathShort", "PathLong", "pathCommon", \
                                     "del_rxn1", "del_rxn2", "totalFluxDiff",\
                                     "solStatus"])

    # Get the reactions
    rxns = df["rxns"]
    new_rxns = []
    for i in rxns:
        temp = []
        for j in i:
            temp.append(str(j[0][0]))
        new_rxns.append(temp)
        
    df["rxns"] = new_rxns
    # Check for empty 
    a = np.where(df["rxns"].apply(len) == 0)
    indices = a[0]
    df.drop(index=indices, inplace=True)

    # Get reaction 1
    del_rxn1 = df["del_rxn1"]
    new_del_rxn1 = []
    for i in del_rxn1:
        new_del_rxn1.append(i[0][0][0])
    df["del_rxn1"] = new_del_rxn1

    # Get reaction 2
    del_rxn2 = df["del_rxn2"]
    new_del_rxn2 = []
    for i in del_rxn2:
        new_del_rxn2.append(i[0][0][0])
    df["del_rxn2"] = new_del_rxn2

    # Get flux difference
    diff_flux = df["diff_flux"]
    new_diff_flux = []
    for i in diff_flux:
        temp = []
        for j in i:
            temp.append(j[0])
        new_diff_flux.append(temp)
        
    df["diff_flux"] = new_diff_flux

    # Get absolute flux difference
    abs_diff_flux = df["abs_diff_flux"]
    new_abs_diff_flux = []
    for i in abs_diff_flux:
        temp = []
        for j in i:
            temp.append(j[0])
        new_abs_diff_flux.append(temp)
        
    df["abs_diff_flux"] = new_abs_diff_flux

    # Get PathShort
    PathShort = df["PathShort"]
    new_PathShort = []
    for i in PathShort:
        temp = []
        for j in i:
            temp.append(str(j[0][0]))
        new_PathShort.append(temp)
        
    df["PathShort"] = new_PathShort
    
    # Get PathLong
    PathLong = df["PathLong"]
    new_PathLong = []
    for i in PathLong:
        temp = []
        for j in i:
            temp.append(str(j[0][0]))
        new_PathLong.append(temp)
        
    df["PathLong"] = new_PathLong
    
    # Get pathCommon
    pathCommon = df["pathCommon"]
    new_pathCommon = []
    for i in pathCommon:
        temp = []
        for j in i:
            temp.append(str(j[0][0]))
        new_pathCommon.append(temp)
        
    df["pathCommon"] = new_pathCommon
    
    # Get totalFluxDiff
    totalFluxDiff = df["totalFluxDiff"]
    new_totalFluxDiff = []
    for i in totalFluxDiff:
        new_totalFluxDiff.append(i[0][0])
        
    df["totalFluxDiff"] = new_totalFluxDiff

    # Get solStatus
    solStatus = df["solStatus"]
    new_solStatus = []
    for i in solStatus:
        new_solStatus.append(i[0][0])
        
    df["solStatus"] = new_solStatus

    sl_size = [len(i) for i in df["rxns"]]
    df["sl_size"] = sl_size

    common_sl_size = [len(i) for i in df["pathCommon"]]
    df["common_sl_size"] = common_sl_size
    df["num_diff"] = get_diff(df["diff_flux"])
    df["net_diff"] = get_net_diff(df["diff_flux"])

    return df

def get_data(model):
    df = pd.read_csv("../results/" + model + '/' + model + '_pFBA.csv')
    df = df.replace(["GLCtex_copy1", "GLCtex_copy2"], "GLCtex")
    rxns = list(df['Rxn_1'].unique())
    rxns.extend(list(df['Rxn_2'].unique()))
    rxns = list(set(rxns))
    types = {}
    
    rxns = [i.split("_copy")[0] for i in rxns]
    
    for i in tqdm(rxns, desc=model):
        url = 'http://bigg.ucsd.edu/models/' + model + '/reactions/' + i
        source = requests.get(url).text
        soup = bs(source, 'lxml')

        values = {}
        for data in soup.find_all('div', class_='col-lg-8'):
            all_h4 = [i.text for i in data.find_all('h4')]
            all_p = [i.text for i in data.find_all('p')]
            all_h4.remove('Metabolites:')

            for h4,p in zip(all_h4, all_p):
                values[h4[:-1]] = p
            values[h4[:-1]] = p.split("\n")[1]
        
        if 'Subsystem' not in values.keys():
            print(i, url, values)
        types[i] = values['Subsystem']
    
    rxn1_temp = df["Rxn_1"].str.split("_copy").str[0]
    rxn2_temp = df["Rxn_2"].str.split("_copy").str[0]
    
    subsystems1 = [types[i] for i in rxn1_temp]       
    subsystems2 = [types[i] for i in rxn2_temp]
    df["Subsystem 1"] = subsystems1
    df["Subsystem 2"] = subsystems2
    
    return df

def convert_csv_to_tex(fname, column_headings):
    text_in_tex = "\\begin{table}[ht]\n\\centering\n\\begin{tabular}{l c c c}\n\\hline\n\\hline"
    for i in column_headings[:-1]:
        text_in_tex += "\\textbf{"+i+"} & "
    text_in_tex += "\\etxtbf{"+column_headings[-1]+"}\\\\"

    text_in_tex += "\n\\hline\n\\hline\n"
    
