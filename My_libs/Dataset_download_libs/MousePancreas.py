# Lib.py

import pandas as pd
import numpy as np
from collections import Counter

import anndata
import scanpy


def download_data(data_file_path = '', data_file_name = '', 
                  verbose = True):
    
    print("Welcome to PANCREAS EMBRYOGENESIS dataset!")

    adata_raw=scanpy.read_h5ad(data_file_path + data_file_name)

    mtx_rawcounts = adata_raw.X
    df_meta = adata_raw.obs

    df_meta['cell'] = np.arange(0, len(df_meta), 1)
    df_meta['stage'] = df_meta['stage'].values.astype(str)
    cells = df_meta.cell.values

    if(verbose):
        print(f"Metadata in a dataframe with shape ({len(df_meta)}, {len(df_meta.columns)})")
        print(f"scRNA-seq data in a counts matrix with shape ({mtx_rawcounts.shape})")

    mtx = mtx_rawcounts.tocsr()[cells, :]
    del mtx_rawcounts

    adata_raw.var['gene_name'] = np.array(adata_raw.var.index)
    genes_names = adata_raw.var.gene_name.values

    # ---------------------------------------------------FILTRO CELLULE----------------------------------------------------
    
    if(verbose): print("\nQuality control on cells...")

    N_zeros_cells = np.sum(mtx.getnnz(1) == 0)
    min_occurrence = 1200
    N_sparse_cells = np.sum(mtx.getnnz(1) < min_occurrence)
    cells_cond1 = mtx.getnnz(1) > min_occurrence
    if(verbose):    
        print("In order to follow the quality control of the paper:")
        print(f" - cells with less than {min_occurrence} expressed genes were deleted ({N_sparse_cells})")
    
    MT_ratios = adata_raw.obs.mt_frac.values
    max_fraction = 0.20
    cells_cond2 = MT_ratios < max_fraction
    if(verbose):    print(f" - cells with mitochondrial gene-expression fractions greater than {100*max_fraction}% ({np.sum(cells_cond2==False)}) were deleted")

    cells_cond = (cells_cond1 & cells_cond2)
    
    # -------------------------------------------------------------------------------------------------------
    
    mtx = mtx[cells_cond, :]
    df_meta = df_meta[cells_cond]
    
    # ------------------------------------------------------FILTRO GENI----------------------------------------------------
    
    if(verbose): print("\nGenes selection...")
    
    protCoding_genes = get_protCoding_genes(data_file_path, genes_names)
    genes_cond1 = protCoding_genes
    if(verbose): print(f"Selecting {np.sum(protCoding_genes)} protein-coding genes")

    genes_cond2 = mtx.getnnz(0) > 0
    if(verbose): print("Deleting genes because full of zeros")

    genes_cond = (genes_cond1 & genes_cond2)

    # --------------------------------------------------------------------------------------------------------------------
    
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

    df_meta.reset_index(inplace=True)
    del df_meta['index']
    
    return mtx, df_meta, genes_names






def prepare_data(df_meta, mtx, genes_name, 
                 labels, group, verbose = True):

    df, cells = cells_to_take(group, df_meta, labels, verbose)

    mtx = mtx[cells, :]

    if(verbose): print(f"Sub-sampled data in a csr matrix with shape ({mtx.shape})")

    del df['cell']
    df.insert(0, 'cell', np.arange(0, len(df), 1))
    df.reset_index(inplace=True)
    del df['index']

    return mtx, df, genes_name





def get_protCoding_genes(data_file_path, genes_names):
    protCoding_filename = data_file_path + "PC_mmusculus_gene_ensembl.csv"

    ProtCoding_df = pd.read_csv(protCoding_filename)
    ProtCoding_names = ProtCoding_df.external_gene_name.values
    mask = np.array([True if g in ProtCoding_names else False for g in genes_names])

    return mask





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
    




def cell_in_smaller_stage(df):
    temp = df.stage.values
    celltype_count = Counter(temp)
    celltype_count = dict(sorted(celltype_count.items(), key=lambda item: item[1], reverse=True))
    values = list(celltype_count.values())
    
    return values[-1]



def cell_in_smaller_celltype(df):
    temp = df.celltype.values
    celltype_count = Counter(temp)
    celltype_count = dict(sorted(celltype_count.items(), key=lambda item: item[1], reverse=True))
    values = list(celltype_count.values())
    
    return values[-1]