import pandas as pd

def cal_hUSI(adata):
    mm_l2 = pd.read_csv('/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.csv',index_col=0)
    genes = set(mm_l2.index) & set(adata.var_names)
    try:
        exp = adata[:,list(genes)].X.todense()
    except:
        exp = adata[:,list(genes)].X
    exp = pd.DataFrame(exp,index=adata.obs_names,columns=list(genes))
    score = []
    for row in range(len(exp)):  
        score.append(mm_l2.w[genes].corr(exp.iloc[row],method='spearman'))
    return score