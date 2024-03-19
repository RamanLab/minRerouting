import scipy.io
import numpy as np
from analysis.utils import *
import pandas as pd
import seaborn as sns

model_list = ['e_coli_core', 'iECW_1372', 'iIT341', 'iJO1366', 'iML1515', 'iNJ661', 'iPC815', 'iSSON_1240', 'iYL1228', 'iYS1720', 'STM_v1_0']

norm = "two"
df_main, df_select, psl_df_main = consolidate_results(model_list, norm)
consolidate_results_print_dfs(df_main, norm)