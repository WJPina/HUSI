import pandas as pd

mm_l2 = pd.read_csv('/home/wangjing/wangj/AgingScore/Data/Bulk_TrainModel/mm_l2.csv',index_col=0)

def cal_hUSI(adata,mm_l2):
    genes = set(mm_l2.index) & set(adata.var_names)
    exp = adata[:,list(genes)].X.todense()
    exp = pd.DataFrame(exp,index=adata.obs_names,columns=list(genes))
    score = []
    for row in range(len(exp)):  
        score.append(mm_l2.w[genes].corr(exp.iloc[row],method='spearman'))
    return score