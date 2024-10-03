# Lib.py

import pandas as pd
import numpy as np
from collections import Counter

import anndata
import scanpy



def download_data(data_file_path = '', data_file_name = '', 
                  verbose = True):

    print("Welcome to ZEBRAFISH RAJ dataset!")

    adata_raw=scanpy.read(data_file_path + data_file_name)

    df_meta = adata_raw.obs
    mtx_rawcounts = adata_raw.X
    df_meta['cell'] = np.arange(0, len(df_meta), 1)
    df_meta['stage'] = df_meta['stage'].values.astype(str)
    cells = df_meta.cell.values

    if(verbose):
        print(f"scRNA-seq data in a counts matrix cells x genes with shape ({mtx_rawcounts.shape})")
        print("Gene names stored in adata.var")
        print(f"Metadata about cells stored in adata.obs ({df_meta.columns})")

    mtx = mtx_rawcounts.tocsr()[cells, :]
    del mtx_rawcounts

    genes_names = adata_raw.var.gene_name.values

    # -------------------------------------------Filter on cells-------------------------------------------
    
    if(verbose): print("\nQuality control on cells...")

    # --------------------------------------------------Filter on genes----------------------------------------------
    
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
    
    # ---------------------------------------------------------------------------------------------------------

    if(verbose): print(f"\nNormalization of the counts matrix...")
    adata = anndata.AnnData(mtx)
    scanpy.pp.normalize_total(adata, inplace=True, target_sum=1.)
    mtx = adata.X

    del df_meta['cell']
    df_meta.insert(0, 'cell', np.arange(0, len(df_meta), 1))

    if(verbose): 
        print(f"\nAfter the filtering procedure, scRNA-seq data have shape ({mtx.shape})")
    
    return mtx, df_meta, genes_names





def prepare_data(df_meta, mtx_rawcounts, genes_name, 
                 labels, group, verbose = True):

    df, cells = cells_to_take(group, df_meta, labels, verbose)

    mtx = mtx_rawcounts[cells, :]

    if(verbose): print(f"Sub-sampled data in a matrix with shape ({mtx.shape})")

    del df['cell']
    df.insert(0, 'cell', np.arange(0, len(df), 1))

    return mtx, df, genes_name






def get_protCoding_genes(data_file_path, genes_names):
    protCoding_filename = data_file_path + "PC_drerio_gene_ensembl.csv"

    ProtCoding_df = pd.read_csv(protCoding_filename)
    ProtCoding_names = ProtCoding_df.ensembl_gene_id.values.astype(str)
    ProtCoding_names = np.array([g.lower() for g in ProtCoding_names])
    mask = np.array([True if g.lower() in ProtCoding_names else False for g in genes_names])

    return mask





def cells_to_take(group, df_meta, labels, verbose):

    if(group=="Time"):
        stages_selected = np.array([True if (t in labels) else False for t in df_meta.stage.values])
        df = df_meta[stages_selected]

        n_cells_max = int(cell_in_smaller_stage(df)*0.75)
        n_cells_max = np.min([n_cells_max, 5000])
        print(f"\nSub-sampling so that every stage is equally represented by {n_cells_max} cells")
        df = df.groupby("stage", observed = True).sample(n = n_cells_max)
        cells = df.cell.values

    if(group=="Celltype"):
        celltypes_selected = np.array([True if (ct in labels) else False for ct in df_meta.celltype.values])
        df = df_meta[celltypes_selected]

        n_cells_max = int(cell_in_smaller_celltype(df)*0.75)
        n_cells_max = np.min([n_cells_max, 5000])
        print(f"\nSub-sampling so that every celltype is equally represented by {n_cells_max} cells")
        df = df.groupby("celltype", observed = True).sample(n = n_cells_max)
        cells = df.cell.values
    
    return df, cells






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