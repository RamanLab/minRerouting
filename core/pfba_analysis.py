import os
import pandas as pd
from pathlib import Path
import plotly.express as px
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = 8,6

import warnings
warnings.filterwarnings("ignore")

# two_up =  os.path.abspath(os.path.join(os.getcwd(), "../..", "examples"))
# os.chdir(two_up)

def analyse_pFBA_results(fname, opname):
    print("\n", "="*80, sep="")
    print('Analysing results of', opname, "...", end="")

    df = pd.read_csv(fname)
    grouped = df.groupby(['Rxn_1_Class', 'Rxn_2_Class']).count()
    grouped['Rxn_2'] /= grouped['Rxn_2'].sum()
    grouped = grouped.rename(columns={'Rxn_1':'Count', 'Rxn_2':'Fraction'})
    grouped = grouped.reset_index()
    grouped["Types"] = grouped["Rxn_1_Class"].str.cat(grouped["Rxn_2_Class"],sep=", ")
    
    fig = px.pie(grouped, values='Fraction', names='Types', title='Reaction pair distribution')
    
    grouped.to_csv("../results/"+opname+"/"+opname+"_pfBA_count.csv")
    image_name = "../results/images/"+opname+"_pFBA_count.png"
    fig.write_image(image_name)

    print("Done!")
    print("="*80, sep="")

    image = plt.imread(image_name)
    ax = plt.gca()
    ax.imshow(image)
    title = opname + " pFBA distribution"
    ax.set_title(title)
    ax.axis('off')
    
    return grouped, ax

##########################################################
# Analysis for different models
##########################################################
for model in ['iIT341', 'iML1515', 'iNJ661', 'iPC815', 'iYL1228', 'STM_v1_0', 'e_coli_core']:
    model_results, ax = analyse_pFBA_results("../examples/"+model+"/"+"/"+model+"_pFBA.csv", opname=model)
    print(model_results)
    print("\n", "="*80, sep="")
    plt.show()
