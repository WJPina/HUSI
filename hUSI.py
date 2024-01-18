#######################################################
# This script is used to calculate the hUSI score and
# fit the hUSI score by gassian mixtrue model in python.
#######################################################

import pandas as pd
from tqdm import tqdm
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
import numpy as np

def cal_hUSI(adata):
    mm_l2 = pd.read_csv('./mm_l2.csv',index_col=0)
    genes = set(mm_l2.index) & set(adata.var_names)
    try:
        exp = adata[:,list(genes)].X.todense()
    except:
        exp = adata[:,list(genes)].X
    exp = pd.DataFrame(exp,index=adata.obs_names,columns=list(genes))
    score = []
    for row in tqdm(range(len(exp))):
        score.append(mm_l2.w[genes].corr(exp.iloc[row],method='spearman'))
    return score

def fit_hUSI(score,g=None):
    robjects.r('set.seed(223)')
    mclust = importr('mclust')

    data = list(score)
    data = [np.log2(1+i)/np.log2(1-i) for i in data]
    r_data = FloatVector(data)
    if g is None:
        result = mclust.Mclust(r_data)
    else:
        result = mclust.Mclust(r_data,G=g)  
    age_class = result.rx2('classification')
    age_class = list(age_class)
    return age_class