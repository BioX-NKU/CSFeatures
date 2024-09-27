# CSFeatures: Identification of Cell-Type-Specific Differential Features in Single-Cell and Spatial Omics Data


## Introduction
CSFeatures is a tool designed to identify differentially expressed genes or differential accessibility regions in single-cell and spatial data.

## Installation
It is recommended to create a new virtual environment using conda to run this project.

```
conda create -n CSFeatures python=3.10
conda activate CSFeatures
pip install -r requirements.txt
```

## Quick Start

### 1. Prepare Data
The data to be prepared includes an expression matrix for scRNA-seq or scATAC-seq (with rows as cells and columns as genes) and the corresponding cell types. For spatial data, spatial coordinates must also be provided. 

In Python, you need to provide an AnnData object as input. AnnData is a class designed to store single-cell data.

The provided AnnData format should meet the following requirements:
- The shape of AnnData should be (number of cells, number of genes or regions).
- The AnnData.obs must contain a 'celltype' field to provide cell classifications.
- If using spatial scRNA-seq/scATAC-seq data, the AnnData.obsm should include a 'spatial' field to provide spatial information.

The following is an example of scRNA-seq/scATAC-seq data:

```bash
AnnData object with n_obs × n_vars = 1000 × 2000
    obs: 'celltype'
```
The following is an example of spatial scRNA-seq/scATAC-seq data:

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

This repository provides four datasets in the `example_data` folder: scATAC-seq, scRNA-seq, spatial scATAC-seq, and spatial scRNA-seq data. 

The following sections demonstrate the complete workflow for using this tool to identify differential genes for each of these four datasets.

Note that the scATAC-seq and spatial scATAC-seq data require preprocessing before use. We have provided preprocessing scripts for this purpose.

```python
python atac_processed.py --input atac_data_path --fpeak 0.05 --output ./
```

## Example:

For scRNA-seq data:
```python
import marker_utils 
import pandas as pd

info, adata = marker_utils.getMarkersEI(adata)  # Use default parameters to identify differentially expressed genes.
print(pd.DataFrame(info['Naive_CD4_T']).head(6))  # Print information for the top six differential genes.
save_data(info, output_dir="./output")  # Save the data to the output folder.
```
The output will be:
|   Raw index   |     mean |      var |   mean_hat |   var_hat |   smoothness | gene_name       |      V |       prop_sum |   prop_max |      prop |   other_clusters_local_mean_max |           EI |
|------|----------|-----------|------------|-----------|--------------|------------------|--------|----------------|------------|-----------|-----------------------------------|---------------|
| 11147| 0.850457 | 1.024728  | 0.882790   | 0.216465  |     0.007735 | CCR7            | 93.509598| 1.393134e+15   |   0.265525 | 6.053247e+18 | 1.374805                          | 1.510977e-03  |
| 8454 | 0.276653 | 0.431631  | 0.292217   | 0.048213  |     0.015024 | NELL2           | 5.094197 | 3.178848e+04   |   0.077088 | 5.483430e+06 | 0.588101                          | 7.017379e-06  |
| 2681 | 0.348172 | 0.557063  | 0.343036   | 0.078822  |     0.014723 | FHIT            | 8.233491 | 5.213285e+06   |   0.117773 | 1.120530e+08 | 0.760505                          | 1.205050e-06  |
| 8355 | 0.174513 | 0.320477  | 0.185690   | 0.054763  |     0.017916 | RP11-291B21.2   | 1.699847 | 1.937045e+02   |   0.036403 | 8.533238e+03 | 0.508250                          | 6.396036e-07  |
| 4219 | 0.176416 | 0.296457  | 0.141884   | 0.034865  |     0.018239 | ADTRP           | 1.706381 | 1.145434e+03   |   0.053533 | 1.750278e+04 | 0.487070                          | 2.265048e-07  |
| 7984 | 0.135792 | 0.208369  | 0.128272   | 0.021249  |     0.018304 | PLEKHB1         | 1.007427 | 1.175668e+03   |   0.053533 | 3.603479e+03 | 0.434600                          | 2.791461e-08  |


for scATAC-seq data：
```python
import marker_utils
info, adata = marker_utils.getMarkersEI(adata)  # Use default parameters to identify differential accessibility regions.
print(pd.DataFrame(info['AC']).head(6))  # Print information for the top six differential regions.
save_data(info, output_dir="./output")  # Save the data to the output folder.
```
The output will be:
| Raw index  |   mean   |    var    |  mean_hat  |   var_hat   | smoothness | regions |     V     |   prop_sum   |  prop_max  |      prop      | other_clusters_local_mean_max |    EI    |
|--------|----------|-----------|------------|-------------|------------|-----------|-----------|--------------|------------|----------------|-------------------------------|----------|
| 83745  | 0.350000 | 0.227500  | 0.316087   | 0.026666    | 0.064868   | chr3:50421761-50422261     | 1.888457  | 5.746579     | 0.005581   | 1.586013e+15   | 0.133333                     | 0.079263 |
| 94812  | 0.325000 | 0.219375  | 0.300165   | 0.026877    | 0.067811   | chr4:96719126-96719626     | 1.557634  | 154.545778    | 0.038140   | 1.301879e+14   | 0.400000                     | 0.000175 |
| 104963 | 0.333333 | 0.222222  | 0.312395   | 0.015747    | 0.064387   | chr5:115436434-115436934    | 1.725673  | 893.704854    | 0.047442   | 2.995592e+14   | 0.333333                     | 0.000081 |
| 46912  | 0.291667 | 0.206597  | 0.263222   | 0.011424    | 0.073111   | chr15:30693733-30694233     | 1.163563  | 26.634278     | 0.021359   | 4.644323e+12   | 0.133333                     | 0.000031 |
| 127974 | 0.308333 | 0.213264  | 0.292535   | 0.020187    | 0.068820   | chr8:93917767-93918267    | 1.381423  | 322.080783    | 0.039683   | 2.458932e+13   | 0.266667                     | 0.000015 |
| 52460  | 0.266667 | 0.195556  | 0.276234   | 0.022122    | 0.073794   | chr16:18083366-18083866     | 0.963640  | 7.299099      | 0.007767   | 3.812292e+11   | 0.133333                     | 0.000008 |


For spatial scRNA-seq data:
```python
import marker_utils
info, adata = marker_utils.get_spatial_MarkersEI(adata)  # Use default parameters to identify differentially expressed genes. Note that the 'spatial' field in adata.obsm provides spatial information.
print(pd.DataFrame(info['Cortex_1']).head(6))  # Print information for the top six differential genes.
save_data(info, output_dir="./output")  # Save the data to the output folder.
```
The output will be:
|   Raw index   |     mean |      var |   mean_hat |   var_hat |   smoothness | gene_name   |      V |       prop_sum |   prop_max |      prop |   other_clusters_local_mean_max |           EI |
|------|----------|-----------|------------|-----------|--------------|--------------|--------|----------------|------------|-----------|-----------------------------------|---------------|
| 5144 | 0.767606 | 0.896695  | 0.760068   | 0.725999  |     0.013568 | Ptgfrn      | 43.428662| 1.120684e+19   |   0.329268 | 2.120345e+22 | 1.133333                          | 0.000237      |
| 20105| 0.869718 | 1.380913  | 0.856605   | 1.141745  |     0.013621 | Tcap        | 55.534137| 2.986091e+19   |   0.341463 | 1.491028e+22 | 0.866667                          | 0.000104      |
| 1234 | 0.739437 | 1.030696  | 0.730297   | 0.844052  |     0.015998 | Lamc2       | 34.177448| 3.140146e+17   |   0.303502 | 3.099983e+20 | 1.266667                          | 0.000085      |
| 20472| 0.239437 | 0.294783  | 0.235277   | 0.237290  |     0.036353 | BC006965    | 1.577033 | 5.261375e+04   |   0.081967 | 1.810139e+08 | 0.333333                          | 0.000032      |
| 23466| 0.274648 | 0.333019  | 0.270495   | 0.268728  |     0.033271 | Hspb3       | 2.267151 | 4.074289e+06   |   0.115854 | 1.497047e+09 | 0.333333                          | 0.000005      |
| 9198 | 0.133803 | 0.137026  | 0.131427   | 0.109803  |     0.041256 | Ocm         | 0.433959 | 3.467788e+02   |   0.036885 | 2.250198e+05 | 0.266667                          | 0.000002      |

For spatial scATAC-seq data:
```python
import marker_utils
info, adata = marker_utils.get_spatial_MarkersEI(adata)  # Use default parameters to identify differential regions. Note that the 'spatial' field in adata.obsm provides spatial information.
print(pd.DataFrame(info['Cartilage_3']).head(6))  # Print information for the top six differential accessibility regions.
save_data(info, output_dir="./output")  # Save the data to the output folder.
```
The output will be:
|   Raw index    |    mean   |    var    |  mean_hat  |   var_hat   | smoothness |              regions              |     V     |    prop_sum     |  prop_max  |      prop      | other_clusters_local_mean_max |    EI    |
|-------|-----------|-----------|------------|-------------|------------|-------------------------------------|-----------|-----------------|------------|----------------|-------------------------------|----------|
| 20225 | 0.844071  | 1.789142  | 0.682070   | 0.163087    | 0.050133   | chr2:145881298-145881799            | 14.211410 | 2.028193e+07    | 0.127119   | 4.866924e+12   | 1.384429                     | 0.023072 |
| 15403 | 0.696012  | 1.699827  | 0.524912   | 0.169244    | 0.055302   | chr17:56007029-56007530             | 8.759784  | 4.062194e+04    | 0.076923   | 5.746806e+09   | 1.303596                     | 0.009059 |
| 25884 | 0.762636  | 1.758128  | 0.576857   | 0.112532    | 0.053648   | chr5:33228289-33228790              | 10.841203 | 2.411353e+06    | 0.105932   | 1.672401e+11   | 1.221360                     | 0.005955 |
| 20967 | 0.787216  | 2.091659  | 0.555961   | 0.126177    | 0.053523   | chr2:173377888-173378389            | 11.578457 | 1.314712e+06    | 0.105263   | 3.100156e+10   | 1.236590                     | 0.002130 |
| 20866 | 0.752692  | 1.890911  | 0.561031   | 0.159778    | 0.056592   | chr2:168451683-168452184            | 10.011061 | 1.016632e+06    | 0.105263   | 1.767657e+10   | 0.810277                     | 0.002073 |
| 7217  | 0.788245  | 1.992646  | 0.602321   | 0.150687    | 0.054260   | chr12:16151031-16151532             | 11.450917 | 3.021226e+06    | 0.107692   | 5.437124e+10   | 1.154678                     | 0.001743 |
