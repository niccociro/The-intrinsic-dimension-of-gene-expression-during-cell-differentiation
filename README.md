# The-intrinsic-dimension-of-gene-expression-during-cell-differentiation
In this repository, we supply five notebooks for the reproduction of some of the results exposed in the paper "The intrinsic dimension of gene expression during cell differentiation".

We reproduce the steps of our analysis, from the raw data stored in an AnnData object *adata* to the plot of the intrinsic dimension (ID) trend.

*adata* essentially stores three objects:
1. the **counts** matrix in *adata.X*
2. a dataframe containing metadata about single genes in *adata.var*, specifically their names.
3. a dataframe containing metadata about single cells in *adata.obs*. The necessary attribute consists in: 
    - the **temporal stage** from which each cell has been sampled (*stage* column), in case we want to study the ID trend in time (Panel 1-2).
    - the **cell-type** of each cell (*celltype* column), reconstructed by the authors of the corresponding paper, in case we want to study the ID per cell-type (Panel 3-4).

In the folder 'Data' there are 6 files produced with BioMart, a data mining tool that allows to export data from Ensembl database:
- 3 of them contain a list of the mitochondrial genes for 3 different species. Mitochondrial genes are typically needed to replicate the quality control on cells performed by the authors of the papers.
- 3 of them contain a list of the protein-coding genes for 3 different species. Protein-coding genes are needed because it is our only features selection, as explained in the paper.

-------------------------------------------------------------------------------------------------------------
**ACTION REQUIRED BEFORE RUNNING THE CODE**

In order to carry out our analyses, the original data must be formatted into an AnnData object, that contains both the counts matrix and the metadata in the way that has been previously explained. After that, this object must be stored in a *h5ad* file in the folder 'Data'. In the notebooks, we precisely specify the online repository that gives free access to every original dataset and explain what attribute of the metadata we take in consideration and how we renamed it (if necessary).
In order to facilitate this preliminary step, it is possible to download the data already formatted in the correct way to be analyzed by our code from this [Google Drive repository](https://drive.google.com/drive/folders/1bm69GFaq8lcXRjAtxbgQIi2j_H_cX6bi?usp=drive_link).
