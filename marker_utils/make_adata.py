import scanpy as sc
import pandas as pd
from scipy.sparse import csr_matrix
import numpy as np



def csv_to_adata(data:pd.core.frame.DataFrame,celltype:pd.core.frame.DataFrame):
    '''
    data: [cells,genes]
    celltype: [cells,1]
    '''
    adata=sc.AnnData(data)
    adata.X=csr_matrix(adata.X).astype(dtype=np.float32)
    adata.obs['celltype']=celltype
    return adata


