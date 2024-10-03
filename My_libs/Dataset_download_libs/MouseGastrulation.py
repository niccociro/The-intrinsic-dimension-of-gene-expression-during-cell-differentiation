# Lib.py

import pandas as pd
import numpy as np
import collections

import anndata
import scanpy



def download_data(data_file_path = '', data_file_name = '', 
                  verbose = True):
    
    print("Welcome to MOUSE GASTRULATION dataset!")
    
    adata_raw=scanpy.read_h5ad(data_file_path + data_file_name)

    mtx_rawcounts = adata_raw.X
    df_meta = adata_raw.obs

    df_meta['cell'] = np.arange(0, len(df_meta), 1)
    df_meta['stage'] = df_meta['stage'].values.astype(str)

    if(verbose):
        print(f"scRNA-seq data in a counts matrix cells x genes with shape ({mtx_rawcounts.shape})")
        print("Gene names stored in adata.var")
        print(f"Metadata about cells stored in adata.obs ({df_meta.obs.columns})")

    if(verbose):    print("Dropping cells with nan values in metadata")
    df, cells = clean_df_meta(df_meta, verbose)
    del df_meta

    mtx = mtx_rawcounts.tocsr()[cells, :]
    del mtx_rawcounts

    genes_names = adata_raw.var.gene_name.values

    # ---------------------------------------------------Filter on cells-------------------------------------
    
    if(verbose): print("\nQuality control on cells...")
    
    min_occurrence = 1000
    N_sparse_cells = np.sum(mtx.getnnz(1) < min_occurrence)
    cells_cond1 = mtx.getnnz(1) > min_occurrence
    if(verbose):    
        print("In order to follow the quality control of the paper:")
        print(f" - cells with less than {min_occurrence} expressed genes were deleted ({N_sparse_cells}) deleted")
    
    MT_ratios = MT_fraction(mtx, genes_names, data_file_path)
    max_fraction = 0.0237
    cells_cond2 = MT_ratios < max_fraction
    if(verbose):    print(f" - cells with mitochondrial gene-expression fractions greater than {100*max_fraction}% were deleted ({np.sum(cells_cond2==False)})")

    cells_cond = (cells_cond1 & cells_cond2)
    
    # -------------------------------------------------------------------------------------------------------
    
    mtx = mtx[cells_cond, :]
    df = df[cells_cond]
    
    # ----------------------------------------------------Filter on genes----------------------------------------
    
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
    
    df.reset_index(inplace=True)
    del df['index']

    return mtx, df, genes_names




def prepare_data(df_meta, mtx_rawcounts, genes_name, 
                 labels, group, verbose):

    df, cells = cells_to_take(group, df_meta, labels, verbose)

    mtx = mtx_rawcounts[cells, :]

    if(verbose): print(f"Sub-sampled data in a matrix with shape ({mtx.shape})")

    del df['cell']
    df.insert(0, 'cell', np.arange(0, len(df), 1))
    df.reset_index(inplace=True)
    del df['index']

    return mtx, df, genes_name





def pos_repetead_genes(mtx, gene_names):
    new_gene_names = gene_names
    names_dict = collections.Counter(gene_names)
    names_dict = dict(sorted(names_dict.items(), key=lambda kv: kv[1], reverse=True))

    repeated_genes = np.array([key for key, val in names_dict.items() if val > 1])

    pos_to_remove = []
    pos_to_keep = []
    for rg in repeated_genes:
        pos = np.where(new_gene_names==rg)[0]
        pos_to_remove.append(np.array(pos))
        
        sums = []
        for p in pos:
            sums.append(np.sum(mtx[:, p].toarray()))
        pos_to_keep.append(pos[np.where(sums==np.max(sums))[0]])

    pos_to_remove = np.concatenate(pos_to_remove)
    pos_to_remove = np.array(pos_to_remove, dtype=object)
    pos_to_keep = np.array(pos_to_keep, dtype=object)

    pos_to_remove = np.array([p for p in pos_to_remove if (p in pos_to_keep)==False])
    new_gene_names = new_gene_names[np.array([ind for ind, gene in enumerate(new_gene_names) 
                                              if ((ind in pos_to_remove)==False)])]

    return pos_to_remove




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

        if(verbose): print(f"Selecting cells from the same sample for each stage")
        cells = filter_sample(df_meta, group)
        df = df[df_meta['cell'].isin(cells)]

        n_cells_max = int(cell_in_smaller_stage(df)*0.75)
        n_cells_max = np.min([n_cells_max, 5000])
        if(verbose):    print(f"\nSub-sampling so that every stage is equally represented by {n_cells_max} cells")
        df = df.groupby("stage", observed = True).sample(n = n_cells_max)
        cells = df.cell.values

    if(group=="Celltype"):
        celltypes_selected = np.array([True if (ct in labels) else False for ct in df_meta.celltype.values])
        df = df_meta[celltypes_selected]

#--------------------In case of single cell-types every sample is taken-----------------

        n_cells_max = int(cell_in_smaller_celltype(df)*0.75)
        n_cells_max = np.min([n_cells_max, 5000])
        if(verbose):    print(f"\nSub-sampling so that every celltype is equally represented by {n_cells_max} cells")
        df = df.groupby("celltype", observed = True).sample(n = n_cells_max)
        cells = df.cell.values
    
    return df, cells






def filter_sample(df, group):
    if(group == 'Time'):    labels = np.unique(df.stage.values)
    if(group == 'Celltype'):    labels = np.unique(df.celltype.values)
    cells = []

    for label in labels:
        if(group == 'Time'):    selected_indexes = (df.stage.values == label)
        if(group == 'Celltype'):    selected_indexes = (df.celltype.values == label)
        
        temp_dict = collections.Counter(df[selected_indexes]['sample'].values)
        
        temp_dict = dict(sorted(temp_dict.items(), key=lambda item: item[1], reverse=True))
        major_samples = list(temp_dict.keys())[:1]

        if(group == 'Time'):    selected_indexes = (df.stage.values == label) & (np.isin(df['sample'], major_samples))
        if(group == 'Celltype'):    selected_indexes = (df.celltype.values == label) & (np.isin(df['sample'], major_samples))
        cells.append(df.cell.values[selected_indexes])

    cells = np.concatenate(cells)
    return cells










def clean_df_meta(df_meta, verbose):

    df = df_meta
    
    del df['cell']
    df.insert(0, 'cell', np.arange(0, len(df), 1))
    
    # Drop cells from 'mixed_gastrulation' stage (spurious)
    df.drop(df[df.stage == 'mixed_gastrulation'].index, inplace=True)

    # Clean data
    if(verbose): print("\nCleaning meta data...")
    df.replace({'doublet': 'FALSE'}, False)
    # Drop doublets
    df.drop(df[df.doublet == True].index, inplace=True)
    # Drop cells with nan celltype
    df.drop(df[df.celltype == 'nan'].index, inplace=True)
    df['sample'] = df['sample'].astype('int16')
    
    cells = df.cell.values

    del df['doublet']

    del df['cell']
    df.insert(0, 'cell', np.arange(0, len(df), 1))

    return df, cells





def MT_fraction(mtx, genes_name, data_file_path):
    filename = data_file_path + "MT_mmusculus_gene_ensembl.txt"
    df = pd.read_csv(filename)
    MT_genes = np.unique(np.array(df["Gene name"].values))

    mask = np.array([True if (gname in MT_genes) else False for gname in genes_name])

    mtx0_array = mtx.toarray()
    row_sums = np.sum(mtx0_array, axis=1)
    del mtx0_array

    mtx0_array = mtx[:, mask].toarray()
    MTrow_sums = np.sum(mtx0_array, axis=1)

    # Calculate the sum of MT_genes for each cell and compute the MT_fraction
    cell_MT_fractions = MTrow_sums / row_sums
    
    return cell_MT_fractions



def cell_in_smaller_stage(df):
    temp = df.stage.values
    celltype_count = collections.Counter(temp)
    celltype_count = dict(sorted(celltype_count.items(), key=lambda item: item[1], reverse=True))
    values = list(celltype_count.values())
    
    return values[-1]



def cell_in_smaller_celltype(df):
    temp = df.celltype.values
    celltype_count = collections.Counter(temp)
    celltype_count = dict(sorted(celltype_count.items(), key=lambda item: item[1], reverse=True))
    values = list(celltype_count.values())
    
    return values[-1]