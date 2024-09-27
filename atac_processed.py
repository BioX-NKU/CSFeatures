# %%
import scipy
import numpy as np
from sklearn.feature_extraction.text import TfidfTransformer
import scanpy as sc
import episcanpy as epi
from pathlib import Path
import argparse

# %%
# Perform TF-IDF (count_mat: peak*cell)
def tfidf1(count_mat): 
    nfreqs = 1.0 * count_mat / np.tile(np.sum(count_mat,axis=0), (count_mat.shape[0],1))
    tfidf_mat = np.multiply(nfreqs, np.tile(np.log(1 + 1.0 * count_mat.shape[1] / np.sum(count_mat,axis=1)).reshape(-1,1), (1,count_mat.shape[1])))
    return scipy.sparse.csr_matrix(tfidf_mat)

# Perform Signac TF-IDF (count_mat: peak*cell)
def tfidf2(count_mat): 
    tf_mat = 1.0 * count_mat / np.tile(np.sum(count_mat,axis=0), (count_mat.shape[0],1))
    signac_mat = np.log(1 + np.multiply(1e4*tf_mat,  np.tile((1.0 * count_mat.shape[1] / np.sum(count_mat,axis=1)).reshape(-1,1), (1,count_mat.shape[1]))))
    return scipy.sparse.csr_matrix(signac_mat)

def tfidf3(count_mat): 
    model = TfidfTransformer(smooth_idf=False, norm="l2")
    model = model.fit(np.transpose(count_mat))
    model.idf_ -= 1
    tf_idf = np.transpose(model.transform(np.transpose(count_mat)))
    return scipy.sparse.csr_matrix(tf_idf)

# %%
parser = argparse.ArgumentParser(description="")

parser.add_argument("--input", help="data path",default=False)
parser.add_argument("--fpeak", help="fpeak",default=0.05)
parser.add_argument("--lazy", help="PCA,tSNE,umap", default=True)
parser.add_argument("--output",default="")
# 解析参数
args = parser.parse_args()
datapath=args.input
fpeak=float(args.fpeak)
lazy=args.lazy
output=args.output
Path(output).mkdir(parents=True, exist_ok=True)

# %%
data=sc.read_h5ad(datapath)

epi.pp.binarize(data)  
epi.pp.filter_features(data, min_cells=np.ceil(fpeak*data.shape[0]))  
data.X = tfidf2(data.X.T.toarray()).T  
if lazy:
    epi.pp.lazy(data)  # PCA、tSNE、UMAP
print(data)
output_path=f"{output}{Path(datapath).stem}_processed.{Path(datapath).suffix}"
print(f"The processed data has been output to {output_path}")
data.write(output_path)


# %%
