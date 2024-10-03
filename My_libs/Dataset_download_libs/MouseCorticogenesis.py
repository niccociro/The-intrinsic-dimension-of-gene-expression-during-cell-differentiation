# Lib.py

import pandas as pd
import numpy as np
from collections import Counter

import anndata
import scanpy



def download_data(data_file_path = '', data_file_name = '', 
                  verbose = True):
    
    print("Welcome to MOUSE CORTICOGENESIS dataset!")

    adata_raw = anndata.read_h5ad(data_file_path + data_file_name)
    
    df = adata_raw.obs

    mtx_rawcounts = adata_raw.X
    df_meta = df
    df_meta['cell'] = np.arange(0, len(df_meta), 1)
    df_meta['stage'] = df_meta['stage'].values.astype(str)
    print(f"Deleting doublet cells ({np.sum(df_meta.celltype == 'Doublet')}) and low quality cells ({np.sum(df_meta.celltype == 'Low quality cells')})")
    df_meta.drop(df_meta[df_meta.celltype == 'Doublet'].index, inplace=True)
    df_meta.drop(df_meta[df_meta.scrublet_doublet == 'True'].index, inplace=True) 
    df_meta.drop(df_meta[df_meta.scrublet_doublet == 'TRUE'].index, inplace=True) 
    df_meta.drop(df_meta[df_meta.celltype == 'Low quality cells'].index, inplace=True)
    cells = df_meta.cell.values

    if(verbose):
        print(f"scRNA-seq data in a counts matrix cells x genes with shape ({mtx_rawcounts.shape})")
        print("Gene names stored in adata.var")
        print(f"Metadata about cells stored in adata.obs ({df_meta.columns})")

    mtx = mtx_rawcounts.tocsr()[cells, :]
    del mtx_rawcounts

    genes_names = adata_raw.var.gene_name.values

    # --------------------------------------------Filter on cells-----------------------------------------
    if(verbose): print("\nQuality control on cells...")

    min_occurrence = 500
    N_sparse_cells = np.sum(mtx.getnnz(1) < min_occurrence)
    cells_cond1 = mtx.getnnz(1) > min_occurrence
    if(verbose):    
        print("In order to follow the quality control of the paper:")
        print(f" - cells with less than {min_occurrence} expressed genes were deleted ({N_sparse_cells})")
    
    MT_ratios = df.percent_mito.values.astype(dtype=np.float32)
    max_fraction = 0.075
    cells_cond2 = MT_ratios < max_fraction
    if(verbose):    print(f" - cells with mitochondrial gene-expression fractions greater than {100*max_fraction}% were deleted ({np.sum(cells_cond2==False)})")

    cells_cond = (cells_cond1 & cells_cond2)
    
    # -------------------------------------------------------------------------------------------------------
   
    mtx = mtx[cells_cond, :]
    df = df[cells_cond]
   
    # ---------------------------------------------Filter on genes-------------------------------------------
    
    if(verbose): print("\nGenes selection...")

    protCoding_genes = get_protCoding_genes(data_file_path, genes_names)
    genes_cond1 = protCoding_genes
    if(verbose): print(f"Selecting {np.sum(protCoding_genes)} protein-coding genes")

    genes_cond2 = mtx.getnnz(0) > 0
    if(verbose): print("Deleting genes full of zeros")

    genes_cond = (genes_cond1 & genes_cond2)

    # -------------------------------------------------------------------------------------------------------
    
    genes_names = genes_names[genes_cond]
    mtx = mtx[:, genes_cond]
    
    # -------------------------------------------------------------------------------------------------------

    if(verbose): print(f"\nNormalization of the counts matrix...")
    adata = anndata.AnnData(mtx)
    scanpy.pp.normalize_total(adata, inplace=True, target_sum=1.)
    mtx = adata.X

    del df['cell']
    df.insert(0, 'cell', np.arange(0, len(df), 1))

    if(verbose): 
        print(f"\nAfter the filtering procedure, scRNA-seq data have shape ({mtx.shape})")
    
    return mtx, df, genes_names






def prepare_data(df_meta, mtx_rawcounts, genes_name, 
                 labels, group, verbose = True):

    df, cells = cells_to_take(group, df_meta, labels, verbose)

    mtx = mtx_rawcounts[cells, :]

    if(verbose): print(f"Sub-sampled data in a matrix with shape ({mtx.shape})")

    del df['cell']
    df.insert(0, 'cell', np.arange(0, len(df), 1))

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