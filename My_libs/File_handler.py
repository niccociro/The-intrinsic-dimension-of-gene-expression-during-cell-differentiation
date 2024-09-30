# File_handler.py

import os
import csv
import numpy as np
import json




def generate_covEigFileName(group, n_subsamplings):
    
    file_name = f"{group}"
    file_name = file_name + f"_nsubsamplings{n_subsamplings}"
    
    return file_name


def handle_numpy_array(obj):
  if isinstance(obj, np.ndarray):
    return obj.tolist()
  return obj


def write_to_JSON(my_dict, filename):
    with open(filename, "w") as f:
        json.dump(my_dict, f, default=handle_numpy_array)


def handle_decoded_object(obj):
  if isinstance(obj, list) and all(isinstance(item, float) for item in obj):
    return np.array(obj)
  return obj


def read_from_JSON(filename):
    with open(filename, "r") as f:
        data = json.load(f, object_hook=handle_decoded_object)
    
    return data
#----------------------------------------------------------------------------------------------------