# The-intrinsic-dimension-of-gene-expression-during-cell-differentiation
Notebooks for the reproduction of some results exposed in the paper "The intrinsic dimension of gene expression during cell differentiation".

We reproduce the steps of our analysis, from the raw data to the plot of the intrinsic dimension (ID) trend.
It is important to note that we do not report here the preliminary step consisting in arranging the data into an *anndata* object *adata*, stored in a *h5ad* file. 

*adata* essentially stores two objects:
- the **counts** matrix in *adata.X*
- a dataframe containing metadata about single cells in *adata.obs*. For the following results, the necessary piece of information consists in:
    - **temporal stage** from which each cell has been sampled (*stage* column), in case we want to study the ID trend in time.
    - **cell-type** of each cell (*celltype* column), reconstructed by the authors of the corresponding paper, in case we study the ID per cell-type.

For transparency, for each dataset, we specify the online repository that gives free access to the original data and explain what metadata we take in consideration for our analysis.
