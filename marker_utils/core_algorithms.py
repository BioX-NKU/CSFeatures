import episcanpy as epi
from sklearn.neighbors import kneighbors_graph
import numpy as np
from collections import Counter
from sklearn.decomposition import PCA
import pandas as pd
from sklearn.metrics.pairwise import euclidean_distances
import statsmodels.api as sm
from collections import defaultdict
import numpy as np
from scipy.sparse import csr_matrix
from scipy.stats import norm
import scipy
from scipy.optimize import minimize
from scipy.stats import nbinom
from tqdm import tqdm
from joblib import Parallel, delayed
import joblib
import random
import anndata
import numba
from numba import njit,prange,float32
import os
from scipy.stats import poisson
from pathlib import Path
import sys
import math
import scanpy as sc
from copy import copy


def calculate_mean_and_var(adata,cluster_count_matrix:np.ndarray,
    cluster_processed_matrix:np.ndarray,
    gene_name:np.ndarray,
    cluster:str,
    clusters:np.ndarray,
    count_matrix:np.ndarray,
    processed_matrix:np.ndarray,
    ans:dict):
    ans[cluster]["mean"]=np.array(cluster_count_matrix.mean(axis=0)).reshape(-1)

    tmp=copy(cluster_count_matrix)
    tmp.data=tmp.data ** 2
    ans[cluster]["var"]= np.array(tmp.mean(axis=0)) - ans[cluster]["mean"] ** 2
    ans[cluster]["var"]= ans[cluster]["var"].reshape(-1)
    ans[cluster]["mean_hat"]=np.array(cluster_processed_matrix.mean(axis=0)).reshape(-1)

    tmp=copy(cluster_processed_matrix)
    tmp.data=tmp.data ** 2
    ans[cluster]["var_hat"]= np.array(tmp.mean(axis=0)) - ans[cluster]["mean_hat"] ** 2
    ans[cluster]["var_hat"]= ans[cluster]["var_hat"].reshape(-1)

def expression_matrix_correction(adata:anndata._core.anndata.AnnData):
    D=copy(adata.obsp['distances'])
    D.data=1/(D.data+1)
    D=D+scipy.sparse.eye(adata.shape[0])
    D=D/D.sum(axis=1)
    adata.obsm['correctX']=D.dot(adata.X)


def calculate_zero_mu_max(adata,cluster_count_matrix:np.ndarray,
    cluster_processed_matrix:np.ndarray,
    gene_name:np.ndarray,
    cluster:str,
    clusters:np.ndarray,
    count_matrix:np.ndarray,
    processed_matrix:np.ndarray,
    ans:dict):
    '''
    Calculate the maximum mean of the other clusters.
    '''
    mu_max=np.zeros_like(ans[cluster['mean']])
    for index in range(len(ans[cluster]['mean'])):
        other_cluster_max=-1e10
        for j in ans.keys():
            if j==cluster: continue
            other_cluster_max=max(other_cluster_max,ans[j]['mean_hat'][index])
        mu_max[index]=other_cluster_max
    ans[cluster]['mu_max']=mu_max

def calculate_V(adata,cluster_count_matrix:np.ndarray,
        cluster_processed_matrix:np.ndarray,
        gene_name:np.ndarray,
        cluster:str,
        clusters:np.ndarray,
        count_matrix:np.ndarray,
        processed_matrix:np.ndarray,
        ans:dict):
    '''
    Define a method to compute the metrics for a cluster.
    After calculating the metrics, record the metrics using ans[cluster]["metric name"] = metrics array.
    Input.
    cluster_processed_matrix: matrix corresponding to the current cluster. Shape (number of cells, genes).
    cluster_processed_matrix: the processed matrix corresponding to the current cluster. Shape (number of cells, genes).
    gene_name:name of the gene, length is an np array of the number of genes.
    cluster:name of current cluster.
    clusters:list of clusters, length is an np array of cells.
    count_matrix: matrix of expression values, shape (number of cells, genes).
    processed_matrix: matrix of processed values, shape (number of cells, genes).
    ans:Dictionary to record the computed metrics.
    '''
    ans[cluster]["gene_name"]=gene_name
    var=ans[cluster]["var"]
    mean_hat=ans[cluster]["mean_hat"]**2
    # mean_hat=mean_hat/mean_hat.max()
    mean=ans[cluster]["mean"]**2
    # V=mean_hat/((ans[cluster]["smoothness"]))
    V=mean/((ans[cluster]["smoothness"]))
    ans[cluster]["V"]=V


def calculate_local_mean_max(adata,cluster_count_matrix:np.ndarray,
    cluster_processed_matrix:np.ndarray,
    gene_name:np.ndarray,
    cluster:str,
    clusters:np.ndarray,
    count_matrix:np.ndarray,
    processed_matrix:np.ndarray,
    ans:dict):
    '''
    Computes the maximum local mean in the other clusters.
    Will add a local_mean_max to each cluster in the ans dictionary.
    '''
    other_clusters=(clusters==cluster)==False
    X_pca=adata.obsm["X_pca"]
    other_clusters_pca=X_pca[other_clusters] # other cells * genes
    other_clusters_matrix=count_matrix[other_clusters]
    tmp_anndata=sc.AnnData(other_clusters_matrix)
    tmp_anndata.obsm['X_pca']=other_clusters_pca
    epi.pp.neighbors(tmp_anndata,metric='euclidean',n_neighbors=15,use_rep='X_pca')
    knn=tmp_anndata.obsp['distances']+scipy.sparse.eye(tmp_anndata.shape[0])
    knn.data=np.ones_like(knn.data)
    local_mean_matrix=knn.dot(other_clusters_matrix)/knn[0].sum()
    local_mean_max=local_mean_matrix.max(axis=0)
    ans[cluster]["other_clusters_local_mean_max"]=local_mean_max.toarray().reshape(-1)


def calculate_prop(adata,cluster_count_matrix:np.ndarray,
        cluster_processed_matrix:np.ndarray,
        gene_name:np.ndarray,
        cluster:str,
        clusters:np.ndarray,
        count_matrix:np.ndarray,
        processed_matrix:np.ndarray,
        ans:dict,**kwargs):
    '''
    Calculate the percentage of cells with non-zero expression values in the current cluster.
    Input: original expression matrix cluster_count_matrix
    Output: proportion of non-zero clusters other than the current cluster (square and add one to increase their influence weights)
    '''
    other_clusters=list(Counter(clusters).keys())
    other_clusters.remove(cluster)
    prop=np.zeros(shape=[len(other_clusters),len(gene_name)])
    for k,i in enumerate(other_clusters):
        tmp_cluster=clusters==i
        tmp_cluster_matrix=count_matrix[tmp_cluster]
        tmp_cluster_matrix=np.array((tmp_cluster_matrix>0).mean(axis=0)).reshape(-1)
        prop[k,:]=tmp_cluster_matrix
    # prop_sum = np.sum(np.exp(2*prop*100), axis=0)
    prop_sum =np.sum((np.e+1)**(prop*100),axis=0)
    # prop_sum = prop.sum(axis=0)
    # prop_sum = (((prop_sum*100)**0.5)+1)
    ans[cluster]["prop_sum"]=prop_sum
    # ans[cluster]["prop_sum"]=(prop_sum**2)+1
    ans[cluster]["prop_max"]=prop.max(axis=0)
    prop = (cluster_count_matrix>0).mean(axis=0)
    # prop=(prop*100)**2
    prop=np.exp(prop*100)
    ans[cluster]["prop"]=np.array(prop).reshape(-1)


@njit(parallel=True)
def __smoothness(gene_count,cell2cell,similarity_matrix):
    cell_dim=cell2cell.shape[0]
    ans=0
    for j in prange(cell_dim):
        for k in prange(similarity_matrix.shape[1]):
            cell1_index=j
            cell2_index=cell2cell[j,k]
            cell1_count=gene_count[cell1_index]
            cell2_count=gene_count[cell2_index]
            similarity_1and2=similarity_matrix[j,k]
            if ((cell1_count==0 and cell2_count==0)): 
                ans+=1
            else:
                ans+=similarity_1and2*(np.abs(cell1_count-cell2_count)) / (cell1_count+cell2_count+1)
    return ans

from joblib import Parallel, delayed
import numpy as np

def _smoothness(cluster_count_matrix, similarity_matrix):
    '''
    cluster_count_matrix: expression matrix, shape (cells, genes).
    similarity_matrix: similarity matrix, shape (cell, cell).
    '''
    cell_dim = cluster_count_matrix.shape[0]
    gene_dim = cluster_count_matrix.shape[1]

    def process_gene(i):
        gene_count = cluster_count_matrix[:, i].toarray().reshape(-1)
        neighnors = np.max((similarity_matrix > 0).sum(axis=1))
        cell2cell = np.zeros([cell_dim, neighnors], dtype=np.int32)
        for j in range(cell_dim):
            start = similarity_matrix.indptr[j]
            end = similarity_matrix.indptr[j + 1]
            cell2cell[j][:end - start] = similarity_matrix.indices[start:end]
        tmp_similarity_matrix = np.array([similarity_matrix[i, cell2cell[i]].toarray().reshape(-1) for i in range(cell_dim)])
        return __smoothness(gene_count, cell2cell, tmp_similarity_matrix)

    ans = Parallel(n_jobs=16)(delayed(process_gene)(i) for i in range(gene_dim))
    return np.array(ans)

import scipy.sparse

def calculate_smoothness(adata,cluster_count_matrix:np.ndarray,
        cluster_processed_matrix:np.ndarray,
        gene_name:np.ndarray,
        cluster:str,
        clusters:np.ndarray,
        count_matrix:np.ndarray,
        processed_matrix:np.ndarray,
        ans:dict):
    '''
    Calculate the smoothness of each cluster.
    '''
    X_pca=adata.obsm["X_pca"]
    clusters_pca=X_pca[clusters==cluster] # other cells * genes
    tmp_anndata=sc.AnnData(cluster_processed_matrix)
    tmp_anndata.obsm['X_pca']=clusters_pca
    epi.pp.neighbors(tmp_anndata,metric='euclidean',n_neighbors=15,use_rep='X_pca')
    similarity_matrix=tmp_anndata.obsp['distances']
    similarity_matrix.data=1/(similarity_matrix.data+1)
    similarity_matrix= similarity_matrix+scipy.sparse.eye(similarity_matrix.shape[0])
    ans[cluster]["smoothness"]=_smoothness(cluster_count_matrix,similarity_matrix) / similarity_matrix.shape[0]**2
    ans[cluster]["smoothness"]=ans[cluster]["smoothness"]


def calculate_EI(adata,cluster_count_matrix:np.ndarray,
        cluster_processed_matrix:np.ndarray,
        gene_name:np.ndarray,
        cluster:str,
        clusters:np.ndarray,
        count_matrix:np.ndarray,
        processed_matrix:np.ndarray,
        ans:dict,
        **kwargs):
    '''
    Define a method to compute the metrics for a cluster.
    After calculating the metrics, record the metrics using ans[cluster]["metric name"] = metrics array.
    Input:
    cluster_count_matrix: matrix corresponding to the current cluster. Shape (number of cells, genes).
    cluster_processed_matrix: matrix corresponding to current cluster. Shape (number of cells, genes).
    gene_name:name of the gene, length is an np array of the number of genes.
    cluster:name of current cluster.
    clusters:list of clusters, length is an np array of cells.
    count_matrix: matrix of expression values, shape (number of cells, genes).
    processed_matrix: matrix of processed values, shape (number of cells, genes).
    ans:Dictionary to record the computed metrics.
    '''
    # EI=np.array(ans[cluster]['V'])/ ((ans[cluster]["other_clusters_local_mean_max"])**2+1)
    # EI=np.array(ans[cluster]['prop'])/ ((ans[cluster]["prop_sum"]))
    # EI/= ans[cluster]['prop_sum'] 
    # EI*= ans[cluster]['prop']
    
    EI1=np.array(ans[cluster]['V'])/ ((ans[cluster]["other_clusters_local_mean_max"])**2+1)
    EI1=EI1/EI1.max()
    EI2=np.array(ans[cluster]['prop'])/ ((ans[cluster]["prop_sum"]))
    EI2=EI2/EI2.max()
    
    EI=EI1*EI2
    
    # tmp=ans[cluster]["prop_max"]>0.3
    # tmp=tmp.astype(np.float32)
    # tmp[tmp==True]=1e9
    # tmp[tmp==False]=1
    # EI/=tmp 

    ans[cluster]['EI']=EI


import time
is_calculate_time=True
if is_calculate_time:
    time_info_dict=defaultdict(int)
def calculate_gene_info(adata,methods_list=[
    calculate_mean_and_var,
    calculate_smoothness,
    calculate_V,
    calculate_prop,
    calculate_local_mean_max,
    calculate_EI
],celltype_key="celltype"):
    ans=defaultdict(dict)
    genes=np.array(adata.var.index.tolist())
    barcodes=np.array(adata.obs.index.tolist())
    clusters=np.array(adata.obs[celltype_key].tolist())
    count_matrix=adata.X
    processed_matrix=adata.obsm['correctX']
    tqdm0=tqdm(Counter(clusters).keys())
    for cluster in tqdm0:
        for k,method in enumerate(methods_list):
            tqdm0.set_description_str(f"cluster:{cluster}running{method.__name__}({k+1}/{len(methods_list)})")
            if is_calculate_time:
                start_time = time.time()
            method(adata,count_matrix[clusters==cluster],
                processed_matrix[clusters==cluster],
                genes,
                cluster,
                clusters,
                count_matrix,
                processed_matrix,
                ans)
            if is_calculate_time:    
                end_time = time.time()
                time_info_dict[method.__name__]+=end_time-start_time
    if is_calculate_time:
        for i in time_info_dict.keys():
            print(f"{time_info_dict[i]:.2f} seconds used for {i}.")
    return ans


def getMarkersEI(adata:str,n_comps=50,n_neighbors=30,metric='euclidean',method_list:list=[
    calculate_mean_and_var,
    calculate_smoothness,
    calculate_V,
    calculate_prop,
    calculate_local_mean_max,
    calculate_EI
]):
    '''
    input_file:csv file, first column is gene name, the rest of the columns are cells.
    clusters_file:csv file,the first column is cell name, the second column is cell type.
    '''
    # Loading data
    print("------Loading data------")
    adata=sc.read_h5ad(adata)
    print("------Data transposition------")
    print(f"{adata.X.shape[0]} cells, {adata.X.shape[1]} genes.")
    print("------Step1: Calculate PCA------")
    sc.pp.pca(adata,n_comps=n_comps)
    print("------Step2: Calculate the similarity matrix------")
    sc.pp.neighbors(adata,n_neighbors=n_neighbors,use_rep='X_pca',metric=metric)
    print("------Step3: Adjust the expression value------")
    expression_matrix_correction(adata)
    print("------Step4: Calculate EI------")
    info=calculate_gene_info(adata,method_list)
    print("------Congratulations, success!------") 
    for i in info.keys():
        adata.obs[f"{i}_EI"]=info[i]['EI']
    return info,adata


def get_spatial_MarkersEI(adata,n_comps=50,
                        n_neighbors=30,metric='euclidean',
                        spatial_key="spatial",method_list=[
    calculate_mean_and_var,
    calculate_smoothness,
    calculate_V,
    calculate_prop,
    calculate_local_mean_max,
    calculate_EI
]):
    # Loading data
    print("------Loading data------")
    adata=sc.read_h5ad(adata)
    print("------Data transposition------")
    print(f"{adata.X.shape[0]} cells, {adata.X.shape[1]} genes.")
    print("------Step1: Calculate PCA------")
    sc.pp.pca(adata,n_comps=n_comps)
    print("------Step2: Calculate the similarity matrix------")
    sc.pp.neighbors(adata,n_neighbors=n_neighbors,use_rep=spatial_key,metric=metric)
    print("------Step3: Adjust the expression value------")
    expression_matrix_correction(adata)
    print("------Step4: Calculate EI------")
    info=calculate_gene_info(adata,method_list)
    print("------Congratulations, success!------") 
    for i in info.keys():
        adata.var[f"{i}_EI"]=info[i]['EI']
    return info,adata