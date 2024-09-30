import numpy as np
from collections import Counter

from scanpy.pp import highly_variable_genes as hv_genes
import anndata
import scanpy

import gc
import ast



def download_data(data_file_folder = '', verbose = True):
    
    print("Welcome to HYDRA dataset!")

    adata_raw=scanpy.read_h5ad(data_file_folder+'hydra_transcriptome.h5ad')

    df_meta = adata_raw.obs
    df_meta = df_meta.rename(columns={'Cluster': 'celltype'})
    mtx = adata_raw.X
    df_meta['cell'] = np.arange(0, len(df_meta), 1)

    if(verbose):
        print(f"Metadata in a dataframe with shape ({len(df_meta)}, {len(df_meta.columns)})")
        print(f"scRNA-seq data in a counts matrix with shape ({mtx.shape})")

    genes_names = np.array(adata_raw.var.index)

    # -------------------------------------------------FILTRO CELLULE-------------------------------------------
    
    if(verbose): print("\nQuality control on cells...")

    min_occurrence = 200
    max_occurrence = 8000
    cells_cond1 = (mtx.getnnz(1) > min_occurrence) & (mtx.getnnz(1) < max_occurrence)
    N_sparse_cells = np.sum(cells_cond1==False)
    if(verbose):    
        print("In order to follow the quality control of the paper:")
        print(f" - cells with less than {min_occurrence}, or more than {max_occurrence} expressed genes were deleted ({N_sparse_cells}")
    
    cells_size = np.asarray(mtx.sum(axis=1)).flatten()
    min_size = 400
    max_size = 70000
    cells_cond2 = (cells_size > min_size) & (cells_size < max_size)
    if(verbose):    print(f" - cells with size greater than {min_size} and smaller than {max_size}. {np.sum(cells_cond2==False)} deleted.")
    
    MT_ratios = MT_fraction(mtx, genes_names)
    max_fraction = 0.05
    cells_cond3 = MT_ratios < max_fraction
    if(verbose):    print(f" - cells with mitochondrial gene-expression fractions greater than {100*max_fraction}% ({np.sum(cells_cond3==False)}) were deleted")

    cells_cond = (cells_cond1 & cells_cond2 & cells_cond3)

    # -------------------------------------------------------------------------------------------------------
    
    mtx = mtx[cells_cond, :]
    df_meta = df_meta[cells_cond]
    
    # --------------------------------------------------FILTRO GENI----------------------------------------------
    
    if(verbose): print("\nGenes selection...")

    if(verbose): print(f"Cannot select protein-coding genes because they are not available for Hydra in Ensembl database")

    genes_cond = mtx.getnnz(0) > 0
    if(verbose): print("Deleting genes because full of zeros")
    
    # -------------------------------------------------------------------------------------------------------
    
    genes_names = genes_names[genes_cond]
    mtx = mtx[:, genes_cond]
    
    # -------------------------------------------------------------------------------------------------------
    
    if(verbose): print(f"\nNormalization of the counts matrix...")
    adata = anndata.AnnData(mtx)
    scanpy.pp.normalize_total(adata, inplace=True, target_sum=1.)
    mtx = adata.X

    del df_meta['cell']
    df_meta.insert(0, 'cell', np.arange(0, len(df_meta), 1))

    if(verbose): 
        print(f"\nscRNA-seq data in csr matrix with shape ({mtx.shape})")
        print(f"Metadata in a dataframe with columns {list(df_meta.columns)}")

    return mtx, df_meta, genes_names





def prepare_data(df_meta, mtx, genes_name, 
                 labels, group, verbose = True):

    df, cells = cells_to_take(group, df_meta, labels, verbose)

    mtx = mtx[cells, :]

    if(verbose): print(f"Sub-sampled data in a csr matrix with shape ({mtx.shape})")

    del df['cell']
    df.insert(0, 'cell', np.arange(0, len(df), 1))
    df.reset_index(inplace=True)

    return mtx, df, genes_name




def cells_to_take(group, df_meta, labels, verbose):

    if(group=="Time"):
        stages_selected = np.array([True if (t in labels) else False for t in df_meta.stage.values])
        df = df_meta[stages_selected]

        n_cells_max = int(cell_in_smaller_stage(df)*0.75)
        n_cells_max = np.min([n_cells_max, 5000])
        if(verbose):    print(f"\nSub-sampling so that every stage is equally represented by {n_cells_max} cells")
        df = df.groupby("stage", observed = True).sample(n = n_cells_max)
        cells = df.cell.values

    if(group=="Celltype"):
        celltypes_selected = np.array([True if (ct in labels) else False for ct in df_meta.celltype.values])
        df = df_meta[celltypes_selected]

        n_cells_max = int(cell_in_smaller_celltype(df)*0.75)
        n_cells_max = np.min([n_cells_max, 5000])
        if(verbose):    print(f"\nSub-sampling so that every celltype is equally represented by {n_cells_max} cells")
        df = df.groupby("celltype", observed = True).sample(n = n_cells_max)
        cells = df.cell.values
    
    return df, cells
    




def MT_fraction(mtx, genes_name):
    MT_genes_position = []

    for ng, gene in enumerate(genes_name):
        conds = []
        conds.append(gene[:2].lower() == 'mt')
        conds.append(gene[:4].lower() == 'mrpl') #mitochondrial ribosomal proteins
        conds.append(gene[:4].lower() == 'mars') #mitochondrial aminoacyl-tRNA synthetases
        conds.append(gene[:4].lower() == 'mrps') #mitochondrial transcription termination factors
        conds.append(gene[:5].lower() == 'mterf') #mitochondrial ribosomal proteins
        if any(conds):   MT_genes_position.append(ng)
    MT_genes_position = np.array(MT_genes_position)

    # Convert the sparse matrix to a numpy array
    mtx0_array = mtx.toarray()

    # Calculate the sum of each row in the matrix
    row_sums = np.sum(mtx0_array, axis=1)

    # Calculate the sum of MT_genes for each cell and compute the MT_fraction
    #print(MT_genes_position)
    cell_MT_fractions = np.sum(mtx0_array[:, MT_genes_position], axis=1) / row_sums
    
    return cell_MT_fractions




def cell_in_smaller_stage(df):
    temp = df.stage.values
    stage_count = Counter(temp)
    stage_count = dict(sorted(stage_count.items(), key=lambda item: item[1], reverse=True))
    values = list(stage_count.values())
    
    return values[-1]




def cell_in_smaller_celltype(df):
    temp = df.celltype.values
    celltype_count = Counter(temp)
    celltype_count = dict(sorted(celltype_count.items(), key=lambda item: item[1], reverse=True))
    values = list(celltype_count.values())
    
    return values[-1]