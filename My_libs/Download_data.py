import sys

def download_dataset(dataset, data_file_path, data_file_name, verbose):
        
    if(dataset=='ZebraEmbryo_Wagner'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import ZebraEmbryo_Wagner as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.download_data(data_file_path = data_file_path, 
                                                            data_file_name = data_file_name, 
                                                            verbose = verbose)
        

    if(dataset=='ElegansEmbryo'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import ElegansEmbryo as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.download_data(data_file_path = data_file_path, 
                                                            data_file_name = data_file_name, 
                                                            verbose = verbose)
    
    if(dataset=='MouseGastrulation'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import MouseGastrulation as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.download_data(data_file_path = data_file_path, 
                                                            data_file_name = data_file_name, 
                                                            verbose = verbose)
    
    if(dataset=='MousePancreas'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import MousePancreas as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.download_data(data_file_path = data_file_path, 
                                                            data_file_name = data_file_name, 
                                                            verbose = verbose)
    
    if(dataset=='MouseCorticogenesis'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import MouseCorticogenesis as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.download_data(data_file_path = data_file_path, 
                                                            data_file_name = data_file_name, 
                                                            verbose = verbose)
    
    if(dataset=='ZebraNeurogenesis'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import ZebraNeurogenesis as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.download_data(data_file_path = data_file_path, 
                                                            data_file_name = data_file_name, 
                                                            verbose = verbose)
        
    if(dataset=='Hydra'):
        sys.path.insert(1, 'My_libs/Dataset_download_libs/')
        import Hydra as Lib
        import importlib
        importlib.reload(Lib)
        mtx_counts, df_meta, genes_name = Lib.download_data(data_file_path = data_file_path, 
                                                            data_file_name = data_file_name, 
                                                            verbose = verbose)
        
    
    return mtx_counts, df_meta, genes_name

