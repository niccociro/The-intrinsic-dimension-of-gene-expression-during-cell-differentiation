# The-intrinsic-dimension-of-gene-expression-during-cell-differentiation
Notebooks for the reproduction of some results exposed in the paper "The intrinsic dimension of gene expression during cell differentiation".

We reproduce the steps of our analysis, from the raw data stored in an *anndata* object *adata* to the plot of the intrinsic dimension (ID) trend.

*adata* essentially stores two objects:
- the **counts** matrix in *adata.X*
- a dataframe containing metadata about single cells in *adata.obs*. The necessary attribute consists in: 
    - the **temporal stage** from which each cell has been sampled (*stage* column), in case we want to study the ID trend in time (Panel 1-2).
    - the **cell-type** of each cell (*celltype* column), reconstructed by the authors of the corresponding paper, in case we want to study the ID per cell-type (Panel 3-4).

**Warning**
It is important to note that here we do not report the preliminary step consisting in arranging the data into an *anndata* object, stored in a *h5ad* file. Neither we are allowed to redistribute the data in the specific format we need for our analysis.
However, in the notebooks, we specify precisely the online repository that gives free access to every original dataset and explain what metadata we take in consideration for our analysis, in order to give the possibility to repeat autonomously the analysis.
