import os
import pandas as pd
from pathlib import Path
import plotly.express as px

# two_up =  os.path.abspath(os.path.join(os.getcwd(), "../..", "examples"))
# os.chdir(two_up)

def analyse_pFBA_results(fname, opname):
    df = pd.read_csv(fname)
    grouped = df.groupby(['Rxn_1_Class', 'Rxn_2_Class']).count()
    grouped['Rxn_2'] /= grouped['Rxn_2'].sum()
    grouped = grouped.rename(columns={'Rxn_1':'Count', 'Rxn_2':'Fraction'})
    grouped = grouped.reset_index()
    grouped["Types"] = grouped["Rxn_1_Class"].str.cat(grouped["Rxn_2_Class"],sep=", ")
    
    fig = px.pie(grouped, values='Fraction', names='Types', title='Reaction pair distribution')
    
    grouped.to_csv("../../results/"+opname+"/"+opname+"_pfBA_segregation.csv")
    fig.write_image("../../results/images/iML1515_pFBA_count.png")
    
    return grouped, fig

##########################################################
# Analysis for each of the models
##########################################################

iML1515_results, fig = analyse_pFBA_results('../../examples/iML1515/iML1515_pFBA.csv', opname="iML1515")
print(iML1515_results)
fig.show()