# The-intrinsic-dimension-of-gene-expression-during-cell-differentiation
In thi repository, we supply four notebooks for the reproduction of some of the results exposed in the paper "The intrinsic dimension of gene expression during cell differentiation".

We reproduce the steps of our analysis, from the raw data stored in an AnnData object *adata* to the plot of the intrinsic dimension (ID) trend.

*adata* essentially stores two objects:
- the **counts** matrix in *adata.X*
- a dataframe containing metadata about single genes in *adata.var*, specifically their names.
- a dataframe containing metadata about single cells in *adata.obs*. The necessary attribute consists in: 
    - the **temporal stage** from which each cell has been sampled (*stage* column), in case we want to study the ID trend in time (Panel 1-2).
    - the **cell-type** of each cell (*celltype* column), reconstructed by the authors of the corresponding paper, in case we want to study the ID per cell-type (Panel 3-4).

In the folder 'Data' there are 6 files, produced by BioMart, a data mining tool that allows to export data from Ensembl database:
- 3 of them contain a list of the mitochondrial genes for 3 different species. Mitochondrial genes are typically needed to replicate the quality control on cells performed by the authors of the papers.
- 3 of them contain a list of the protein-coding genes for 3 different species. Protein-coding genes are needed because it is our only features selection, as explained in the paper.

**Warning**
In order to carry out the analyses, the original data must be downloaded and formatted into an AnnData object. After that, they have to be stored in a *h5ad* file in the folder 'Data'. In the notebooks, we precisely specify the online repository that gives free access to every original dataset and explain what metadata we take in consideration for our analysis and how we renamed it (if necessary).
