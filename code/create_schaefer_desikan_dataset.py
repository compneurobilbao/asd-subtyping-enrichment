#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script for appending subcortical time series 

@author: javi
"""

from glob import glob
from pathlib import Path
from tqdm import tqdm

output_dir = "/home/javi/Documentos/asd-subtyping-enrichment/data/schaefer_desikan"
Path(output_dir).mkdir(exist_ok=True, parents=True)

ts_filenames_desikan = glob("/home/javi/Documentos/asd-subtyping-enrichment/data/desikan/*")

subcortex_ids = np.array([35,36,37,38,39,40,41,76,77,78,79,80,81,82])
subcortex_ids-=1

for filename in tqdm(ts_filenames_desikan):
    ts_subcortex_desikan = np.loadtxt(filename)[:, subcortex_ids]
    
    filename_schaefer = filename.replace("desikan", "schaefer")
    filename_schaefer = filename_schaefer.replace(".txt", "_schaefer.txt")
    ts_schaefer = np.loadtxt(filename_schaefer)
    
    ts_schaefer_desikan = np.column_stack((ts_schaefer, ts_subcortex_desikan))
    
    filename_schaefer_desikan = filename.replace("desikan", "schaefer_desikan")
    filename_schaefer_desikan = filename_schaefer_desikan.replace(".txt", "_schaefer_desikan.txt")
    np.savetxt(filename_schaefer_desikan, ts_schaefer_desikan)