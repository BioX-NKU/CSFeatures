# CSFeatures: Identification of cell type-specific differential features in single-cell and spatial omics data

### English | [简体中文](README_ZH.md)

## Introduction
CSFeatures is a tool used to identify differential genes or differential peaks in single-cell and spatial data.

## Installation

It is recommended to use conda to create a new virtual environment to run this project.

```
conda create -n CSFeatures python=3.10
conda activate CSFeatures
pip install -r requirements.txt
```

## Quick Start

### 1. Prepare the Data
The data to be prepared includes the expression matrix for RNA or ATAC (with rows as cells and columns as genes) and cell types. For spatial data, spatial coordinate-related data also needs to be provided. Spatial information can be extracted using methods such as STAGATE, SpatialPCA, or SpaGCN.

In Python, you need to provide an Anndata object as input. Anndata is a class designed by anndata for storing single-cell data.

The provided Anndata object should meet the following requirements:

The shape of the Anndata object should be (number of cells, number of genes or peaks).
Anndata.obs should contain a celltype field to provide the cell types.
For spatial RNA/ATAC data, Anndata.obsm should also provide a spatial field to include spatial information.

The following is an example of RNA/ATAC data:
```bash
AnnData object with n_obs × n_vars = 1000 × 2000
    obs: 'celltype'
```
The following is an example of spatial RNA/ATAC data:

```bash
AnnData object with n_obs × n_vars = 1000 × 2000
    obs: 'celltype'
    obsm: 'spatial'
```

### 2. Find Differential Genes/Peaks

This section primarily utilizes the `getMarkersEI` and `get_spatial_MarkersEI` functions.

## `getMarkersEI` Function

This function processes single-cell gene expression data contained in an `AnnData` object and identifies marker genes based on various statistical methods. The function performs the following tasks:
- Dimensionality reduction
- Similarity graph construction
- Expression value adjustment
- Calculation of gene information metrics specified in `method_list`

### Parameters

- **`adata` (str)**:  
  An `AnnData` object containing the single-cell gene expression matrix. Its `.X` attribute should be a 2D array with rows corresponding to cells and columns to genes.

- **`n_comps` (int, optional)**:  
  The number of principal components to calculate in PCA. The default value is 50.

- **`n_neighbors` (int, optional)**:  
  The number of nearest neighbors used to build the k-nearest neighbors graph. The default value is 30.

- **`metric` (str, optional)**:  
  The distance metric used to compute cell similarities. The default is `'euclidean'`.

- **`method_list` (list, optional)**:  
  A list of functions used to compute various gene statistics. The default includes:

  - `calculate_mean_and_var`
  - `calculate_smoothness`
  - `calculate_V`
  - `calculate_prop`
  - `calculate_local_mean_max`
  - `calculate_EI`

  These functions should be defined elsewhere and are used to calculate different metrics for gene expression analysis.

### Return Values

- **`info`**:  
  A data structure (e.g., `DataFrame`) containing gene information metrics calculated based on the methods in `method_list`.

- **`adata`**:  
  The `adata` object with updated EI values.

## `get_spatial_MarkersEI` Function

This function processes single-cell gene expression data containing spatial information and identifies spatial marker genes based on various statistical methods. It performs dimensionality reduction, constructs a similarity matrix, adjusts expression values, and calculates gene information metrics specified in the `method_list` through the following steps.

### Parameters

- **`adata` (AnnData)**:  
  An `AnnData` object containing the single-cell gene expression matrix. Its `.X` attribute should be a 2D array with rows corresponding to cells and columns to genes.

- **`n_comps` (int, optional)**:  
  The number of principal components to calculate in PCA. The default value is 50.

- **`n_neighbors` (int, optional)**:  
  The number of nearest neighbors used to build the k-nearest neighbors graph. The default value is 30.

- **`metric` (str, optional)**:  
  The distance metric used to compute cell similarities. The default is `'euclidean'`.

- **`spatial_key` (str, optional)**:  
  The key used to specify spatial information. The default is `'spatial'`.

- **`method_list` (list, optional)**:  
  A list of functions used to compute various gene statistics. The default includes:

  - `calculate_mean_and_var`
  - `calculate_smoothness`
  - `calculate_V`
  - `calculate_prop`
  - `calculate_local_mean_max`
  - `calculate_EI`

### Return Values

- **`info`**:  
  Data containing gene information metrics calculated based on the methods in `method_list`.

- **`adata`**:  
  The `adata` object with updated EI values.

## Example

If it's RNA/ATAC data, run the following code:
```python
import marker_utils
info,adata=marker_utils.getMarkersEI(adata)
```
If it's spatial RNA/ATAC data, run the following code:
```python
import marker_utils
info,adata=marker_utils.get_spatial_MarkersEI(adata)
```

A sample dataset is provided in the `example_data` folder. You can run the following code to identify differentially expressed genes in the sample data:

```python
import marker_utils
info, adata = marker_utils.get_spatial_MarkersEI("./example_data/pbmc_2700_seurat.h5ad")
print(adata)
```
The output will be:
```bash
AnnData object with n_obs × n_vars = 2700 × 13714
    obs: 'clusters', 'celltype', 'Naive_CD4_T_EI', 'B_EI', 'Memory_CD4_T_EI', 'FCGR3A_Mono_EI', 'NK_EI', 'CD8_T_EI', 'CD14_Mono_EI', 'DC_EI', 'Platelet_EI'
    uns: 'pca', 'neighbors'
    obsm: 'X_pca', 'correctX'
    varm: 'PCs'
    obsp: 'distances', 'connectivities'

```

After running the code, the algorithm will calculate the EI value for each gene in each category. The larger the EI value, the more specific the gene is.

Run the following code to save the results as a CSV file:
```python
save_data(info,output_dir="./output")
```


