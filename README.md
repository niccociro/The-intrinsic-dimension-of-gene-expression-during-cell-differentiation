# The-intrinsic-dimension-of-gene-expression-during-cell-differentiation
Notebooks for the reproduction of some results exposed in the paper "The intrinsic dimension of gene expression during cell differentiation".

We reproduce the steps of our analysis, from the raw data to the plot of the intrinsic dimension (ID) trend.
It is important to note that we do not report here the preliminary step consisting in arranging the data into an *anndata* object *adata*, stored in a *h5ad* file. 

*adata* essentially stores two objects:
- the **counts** matrix in *adata.X*
- a dataframe containing metadata about single cells in *adata.obs*. The necessary attribute consists in:
    - the **temporal stage** from which each cell has been sampled (*stage* column), in case we want to study the ID trend in time (Panel 1-2).
    - the **cell-type** of each cell (*celltype* column), reconstructed by the authors of the corresponding paper, in case we want to study the ID per cell-type (Panel 3-4).

For transparency, in the notebooks, we specify the online repository that gives free access to every original dataset and explain what metadata we take in consideration for our analysis.
