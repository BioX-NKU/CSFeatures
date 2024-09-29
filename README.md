# CSFeatures: Identification of Cell Type-Specific Differential Features in Single-Cell and Spatial Omics Data

## Introduction
CSFeatures is a tool designed to identify cell type-specific differentially expressed genes or differentially accessible regions in single-cell and spatial omics data.

## Installation
It is recommended to create a new virtual environment using conda to run this project.

```
conda create -n CSFeatures python=3.10
conda activate CSFeatures
pip install -r requirements.txt
```

## Quick Start

### 1. Prepare data
The input of CSFeatures consists of a preprocessed gene expression matrix and a vector that contains the corresponding labels for all cells to be investigated, such as cell type labels or clustering labels. For spatial omics data, spatial coordinates must also be provided.

You can also provide an AnnData object as input. The provided AnnData format should meet the following requirements:
- The shape of AnnData should be (number of cells, number of genes or regions).
- The AnnData.obs must contain a 'celltype' field to provide cell classifications.
- If processing spatial RNA-seq or spatial ATAC-seq data, the AnnData.obsm should include a 'spatial' field to provide spatial coordinates.

The following is an example of scRNA-seq/scATAC-seq data:

```bash
AnnData object with n_obs × n_vars = 1000 × 2000
    obs: 'celltype'
```
The following is an example of spatial RNA-seq/ATAC-seq data:

```bash
AnnData object with n_obs × n_vars = 1000 × 2000
    obs: 'celltype'
    obsm: 'spatial'
```

### 2. Find differential features

This section primarily utilizes the `getMarkersEI` and `get_spatial_MarkersEI` functions.

## `getMarkersEI` Function

This function processes gene expression or chromatin accessibility data contained in an AnnData object to identify differential features in single-cell omics data. The function performs the following tasks:
- Dimensionality reduction
- Similarity graph construction
- Calculation of gene information metrics specified in `method_list`

### Parameters

- **`adata` (str)**:  
  An `AnnData` object containing either a single-cell gene expression matrix or chromatin accessibility data. Its `.X` attribute should be a 2D array with rows corresponding to cells and columns to features.

- **`n_comps` (int, optional)**:  
  The number of principal components to calculate in PCA. The default value is 50.

- **`n_neighbors` (int, optional)**:  
  The number of nearest neighbors used to build the k-nearest neighbors graph. The default value is 30.

- **`metric` (str, optional)**:  
  The distance metric used to compute cell similarities. The default is `'euclidean'`.

- **`method_list` (list, optional)**:  
  Function list for calculating various characterization feature activities. The default includes:

  - `calculate_mean_and_var`
  - `calculate_smoothness`
  - `calculate_V`
  - `calculate_prop`
  - `calculate_local_mean_max`
  - `calculate_EI`


### Return Values

- **`info`**:  
 A DataFrame containing feature information metrics calculated based on the methods in `method_list`.

- **`adata`**:  
  The `AnnData` object with updated EI values.

## `get_spatial_MarkersEI` Function

This function processes gene expression or chromatin accessibility data contained in an AnnData object to identify differential features in spatial omics data. The function performs the following tasks:
- Dimensionality reduction
- Similarity graph construction
- Calculation of gene information metrics specified in `method_list`

### Parameters

- **`adata` (AnnData)**:  
   An `AnnData` object containing either a spatial gene expression matrix or chromatin accessibility data. Its `.X` attribute should be a 2D array with rows corresponding to cells and columns to features.
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
 A DataFrame containing feature information metrics calculated based on the methods in `method_list`.

- **`adata`**:  
  The `AnnData` object with updated EI values.

## Tutorial

This repository provides four example datasets on Google Drive: [scRNA-seq](https://drive.google.com/file/d/1LWOnXLHYn8W6GFQ2NTfi84JyY9B4XGOK/view?usp=drive_link),[scATAC-seq](https://drive.google.com/file/d/1mXGWKpOMR4I-mqhyAIQ_UFV6VbHwizdh/view?usp=drive_link),  [spatial scRNA-seq](https://drive.google.com/file/d/1U3_0FIBEcTLzTiAHQG00sMNSLvq7lFtl/view?usp=drive_link) and [spatial scATAC-seq data](https://drive.google.com/file/d/1w7oxnwR_Nma5uTm0yOf4I2O44tGC5Dif/view?usp=drive_link). For detailed workflows on how CSFeatures identifies differentially expressed genes and differentially accessible regions in scRNA-seq, scATAC-seq, spatial RNA-seq, and spatial ATAC-seq data, please refer to the following links:

- [scRNA-seq](./tutorials/scRNA-seq.ipynb)
- [scATAC-seq](./tutorials/scATAC-seq.ipynb)
- [spatial_scRNA-seq](./tutorials/spatial_RNA-seq.ipynb)
- [spatial_scATAC-seq](./tutorials/spatial_ATAC-seq.ipynb)
