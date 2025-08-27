import pandas as pd
import numpy as np

### calculate hUSI scores
def minmax(x):
    return (x-min(x))/(max(x)-min(x))

def cal_hUSI(adata):
    w = pd.read_csv('sc/SenOCLR_features.csv',index_col=0)
    genes = set(w.index) & set(adata.var_names)
    try:
        exp = adata[:,list(genes)].X.todense()
    except:
        exp = adata[:,list(genes)].X
    exp = pd.DataFrame(exp,index=adata.obs_names,columns=list(genes))
    score = []
    for row in range(len(exp)):  
        score.append(w.Weight[list(genes)].corr(exp.iloc[row],method='spearman'))
    score = minmax(score)
    score = pd.Series(score,index=adata.obs_names)
    return score

from concurrent.futures import ThreadPoolExecutor
def calculate_correlation(args):
    row_data, weights = args
    return weights.corr(row_data, method='spearman')

def cal_hUSI_parallel(adata, n_jobs=20):
    # w = model
    w = pd.read_csv('sc/SenOCLR_features.csv',index_col=0)
    genes = list(set(w.index) & set(adata.var_names))
    weights = w.Weight[genes]
    
    try:
        exp = adata[:,genes].X.todense()
    except:
        exp = adata[:,genes].X
    
    exp = pd.DataFrame(exp, index=adata.obs_names, columns=genes)
    row_weight_pairs = [(exp.iloc[i], weights) for i in range(len(exp))]
    
    # Calculate correlations in parallel
    with ThreadPoolExecutor(max_workers=n_jobs) as executor:
        score = list(executor.map(calculate_correlation, row_weight_pairs))
    
    score = minmax(score)
    score = pd.Series(score, index=adata.obs_names)
    
    return score

### SSE
def SSE_hUSI(scores):
    try:
        len(scores.index)
    except:
        message = 'hUSI scores must be a named vector!!!'
        return message
    
    th = None
    min_sse = None
    
    for i in scores:
        low = scores[scores <= i]
        high = scores[scores > i]
        sse = np.sum((low - np.mean(low)) ** 2) + np.sum((high - np.mean(high)) ** 2)
        if th is None and min_sse is None:
            th = i
            min_sse = sse
        elif sse < min_sse:
            th = i
            min_sse = sse
    senclass = np.where(scores <= th, 0, 1)
    senclass=pd.Series(senclass,index=scores.index)
    return senclass

### GMM
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
def GMM_hUSI(scores):
    try:
        len(scores.index)
    except:
        message = 'hUSI scores must be a named vector!!!'
        return message
    
    robjects.r('set.seed(223)')
    mclust = importr('mclust')
    data = scores.tolist()
    r_data = FloatVector(data)
    result = mclust.Mclust(r_data)  
    meanV = result.rx2('parameters').rx2('mean')
    meanV = np.array(meanV)
    meanV = np.sort(meanV)
    scores = scores.sort_values()
    scores_index = scores.index
    senclass = []
    for i in range(len(meanV)):
        class_labels = np.where(scores < meanV[i], i + 1, i + 2)
        if i == len(meanV) - 1:
            senclass.extend(class_labels)
        else:
            senclass.extend(class_labels[class_labels == i + 1])
        scores = scores[class_labels != i + 1]
    senclass = pd.Series(senclass, index=scores_index)
    return senclass

    