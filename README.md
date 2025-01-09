## Using hUSI to evaluate senescence state and reveal novel senescence regulators
![workflow](Fig1.png)
#### We introduced the human universal senescence index (hUSI) into the evaluation of senescence based on transcriptome profiles to identify senescent cells/samples across different cell types, conditions and sequencing platforms. Senescence features were learned from the comprehensive senescence transcriptome profiles (bulk RNA-seq data) utilizing a One-Class Logistic Regression (OCLR) model. The scoring for each sample or cell is determined using the Spearman correlation coefficient. Senescence status can be further classified through Sum of Squares for Error (SSE) or Gaussian Mixture Model (GMM).
## Note
#### We recommend two methods based on different assumptions: (1) If assuming that there are only two cell states, senescence and non-senescence, in the dataset, we used a variance-based method—sum of square error (SSE)—for thresholding. SSE, as the well-performed binary clustering method for one dimension vector, can divide cells into senescence and non-senescence group by minimizing the variance of hUSI of two groups. (2) If assuming there are uncertain senescence states in the dataset, we use Gaussian Mixture Model (GMM)158 to estimate the optimal number of senescence states and the probability that each cell belongs to a specific state. Since GMM is frequently used as optimal thresholding tool in various scenarios, it is helpful to determine cellular senescence heterogeneity.
#### All scriptis for reproducing results in the manuascripts can be found in `HUSI/R` folder. All raw data used for model trianning and validation can be found in `HUSI/Data` folder. All figures in paper can be found in `HUSI/Results` folder.
## Usage
#### Make sure you have already downloaded `HUSI/sc` folder.
#### R
```R
### make sure library packages in classifier.R have been installed
source('classifier.R')
### exp: input normalized gene expression matrix 
library(Seurat)
load('Data/Aarts2017.rdata')
exp = GetAssayData(Aarts2017)
hUSI = cal_hUSI(exp)
### classify senescent group by SSE
SenClass = SSE_hUSI(hUSI)
### classify senescent group by GMM
SenClass = GMM_hUSI(hUSI)
```
#### Python
```python
### make sure import modules in hUSI.py have been installed
from hUSI import cal_hUSI,SSE_hUSI,GMM_hUSI
### adata: input annadata formant with normalized gene expression matrix included as X
import scanpy as sc
adata = sc.read_h5ad('Data/Aarts2017.h5ad')
### make sure library packages in classifier.py and mclust in R have been installed
hUSI = cal_hUSI(adata)
### classify senescent group by SSE
SenClass = SSE_hUSI(hUSI)
### classify senescent group by GMM
SenClass = GMM_hUSI(hUSI)
```
#### package used
`R==4.0.5 Seurat==4.2.1.9001 mclust==6.0.1 dplyr==1.1.4`

`python==3.9.18 pandas==2.1.4 numpy==1.26.3 rpy2==3.5.16`
