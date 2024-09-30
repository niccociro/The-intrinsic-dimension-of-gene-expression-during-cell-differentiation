# DataFiles_Lib.py

import os
import numpy as np
from numpy.linalg import matrix_power
from numpy.random import randint
from dadapy.data import Data
import sys

import File_handler
# Le seguenti 2 righe servono a tenere aggiornata la mia libreria custom se aggiungo cose
import importlib
importlib.reload(File_handler)


#-------------------------------------------------------------------------------------------

def prepare_data(dataset, df, mtx, genes_name, 
                 labels, group, single_organ, 
                 verbose):        
        
    if(dataset=='ZebraEmbryo_Wagner'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import ZebraEmbryo_Wagner as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.prepare_data(df, mtx, genes_name,
                                                                labels, group, verbose)
    
    if(dataset=='ElegansEmbryo'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import ElegansEmbryo as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.prepare_data(df, mtx, genes_name,
                                                            labels, group, verbose)
        
    if(dataset=='MouseGastrulation'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import MouseGastrulation as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.prepare_data(df, mtx, genes_name,
                                                            labels, group, verbose)
    
    if(dataset=='MousePancreas'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import MousePancreas as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.prepare_data(df, mtx, genes_name,
                                                            labels, group, verbose)
    
    if(dataset=='MouseCorticogenesis'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import MouseCorticogenesis as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.prepare_data(df, mtx, genes_name,
                                                            labels, group, verbose)
    
    if(dataset=='ZebraNeurogenesis'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import ZebraNeurogenesis as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.prepare_data(df, mtx, genes_name,
                                                            labels, group, verbose)
        
    if(dataset=='Hydra'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import Hydra as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.prepare_data(df, mtx, genes_name,
                                                            labels, group, verbose)
        
        
    return mtx_counts, df_meta, genes_name







def ID(dataset, group, labels, df, mtx, 
        genes_name, single_organ = 0,
        n_subsamplings = 10, verbose = False):

    impostring = f"Settings:"
    impostring = impostring + f"\nDataset {dataset}"
    impostring = impostring + f"\n{group}-based grouping of cells"
    impostring = impostring + f"\nLabels: {labels}"
    impostring = impostring + f"\n{n_subsamplings} independent sub-samplings"
    impostring = impostring + f"\nCounts matrix dimension={mtx.shape}"
    print(impostring, end = "\n\n")

    overRuns_dict = {}
    seeds = randint(1, 1000000, n_subsamplings)

#-------------------------------------------INIZIO LOOP SUI SUBSAMPLING----------------------------
    
    for n, seed in enumerate(seeds):
    
        print("\nSubsampling number", n+1, ", with seed:", seed, end = ". ")
            
        mtx_rawcounts, df_meta, genes_name0 = prepare_data(dataset, df, mtx, genes_name,
                                                    labels = labels, group = group, 
                                                    single_organ = single_organ,
                                                    verbose = verbose)

#------------------------------------------INIZIO LOOP SU STAGE/CELLTYPE-----------------------------
        
        laio_dict = {}

        for label in labels:
            print(label, end = ' ')

            if(group == 'Time'):    selected_indexes = (df_meta.stage == label)
            elif(group == 'Celltype'):    selected_indexes = (df_meta.celltype == label)
            elif(group == 'Sample'):    selected_indexes = (df_meta.sample_ID == label)
            df_meta_label = df_meta[selected_indexes]
            cells = df_meta_label.cell.values

            mtx1 = mtx_rawcounts[cells, :].toarray()
            
            data = Data(mtx1)
            data.compute_id_2NN(mu_fraction = 0.9)
            laio_dict[label] = data.intrinsic_dim
            print(f"ID = {np.round(data.intrinsic_dim, 1)}", 
                    end = ' - ')

            del mtx1
        
        print('\n')

        overRuns_dict[n] = laio_dict

#----------------------------------------------------------------------------------------------------
        
    file_name = File_handler.generate_covEigFileName(group, n_subsamplings)
 
    folder_measures = "Measures"
    folder_dataset = os.path.join(folder_measures, dataset)

    # Check if "Measures" folder exists
    if not os.path.exists(folder_measures):
        os.makedirs(folder_measures)

    # Check if the folder of the dataset exists
    if not os.path.exists(folder_dataset):
        os.makedirs(folder_dataset)

    measure_file = folder_measures + '/' + dataset + '/' + file_name
    File_handler.write_to_JSON(overRuns_dict, measure_file)

    y_score, y_score_err, (y_shift, y_scaling) = mean_trend(overRuns_dict, 
                                                            labels, n_subsamplings,
                                                            normalize_score=True)

    out_dict = {}
    out_dict['Results from subsamplings'] = overRuns_dict
    out_dict['Mean trend'] = {}
    out_dict['Mean trend']['Labels'] = labels
    out_dict['Mean trend']['Mean'] = y_score
    out_dict['Mean trend']['Std'] = y_score_err
    out_dict['Mean trend']['Normalization params'] = (y_shift, y_scaling)
    
    return out_dict





def mean_trend(overRuns_dict, labels, n_subsamplings, normalize_score):
    y = []
    y_err = []

    for label in labels:
        yr = []

        for subsampling in range(n_subsamplings):
            temp = overRuns_dict[subsampling][label]
            yr.append(temp)

        y.append(np.mean(yr))
        y_err.append(np.std(yr))

    y = np.array(y)
    y_err = np.array(y_err)

    if(normalize_score):
        y_shift = np.min(y-y_err)
        y_scaling = np.max(y+y_err)
        y = y-y_shift
        y = y/(y_scaling-y_shift)
        y_err = y_err/(y_scaling-y_shift)

    return y, y_err, (y_shift, y_scaling)