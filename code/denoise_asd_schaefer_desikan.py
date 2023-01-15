#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 10:33:09 2022

@author: javi
"""

# %% IMPORTS
import numpy as np
import pandas as pd
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import tqdm
from os.path import join as opj
from scipy.io import savemat, loadmat
import os
import pickle
from scipy.stats import pearsonr, ttest_ind, kruskal, levene
from scipy.spatial.distance import squareform, pdist
from statsmodels.stats.multitest import multipletests

# %% PARSE CASE
parser = argparse.ArgumentParser(description='Harmonise and denoise ASD before clustering')
parser.add_argument('--case',
                type=str,
                choices = ["schaefer", "schaefer_desikan"],
                help="Name for the output directory")

opts = parser.parse_args()


# %% READ DATA
data_dir = "/home/javi/Documentos/asd-subtyping-enrichment/data"

# Read ASD subjects
f = open(opj(data_dir, "subjects_subgroups_FD03.pkl"), "rb")
subjects_subgroups = pickle.load(f, encoding='utf-8')
f.close()

asd_subjects = subjects_subgroups['asd_subjs']
td_subjects = subjects_subgroups['tc_subjs']
subjects = np.concatenate((td_subjects, asd_subjects)) #Concatenate subjects to have first TDC and then ASD

# Load demo data
demo = pd.read_excel(opj(data_dir, "demo_final_FD03.xlsx"))

# Compute non interest sources of variability
# Beware, first subjects on the merge to preserve this order
print(np.allclose(subjects, 
                  pd.merge(pd.DataFrame({'SUB_ID': subjects}), demo).SUB_ID))

# Extract different variables (age, sex, site id, group, motion and code ID)
age= pd.merge(pd.DataFrame({'SUB_ID': subjects}), demo).AGE_AT_SCAN.values
sex = pd.merge(pd.DataFrame({'SUB_ID':subjects}), demo ).SEX.values
site = pd.merge(pd.DataFrame({'SUB_ID':subjects}), demo).SITE_ID_MERGED
group = pd.merge(pd.DataFrame({'SUB_ID':subjects}), demo).DX_GROUP
motion = pd.merge(pd.DataFrame({'SUB_ID':subjects}), demo).loc[:, ["FD_mean", "DVARS_mean"]].to_numpy()
code_id = pd.merge(pd.DataFrame({'SUB_ID':subjects}), demo).loc[:, "ID"].to_list()

# Use code ID to extract time series filenames
#ts_filenames = ["ts_" + name + "_schaefer.txt" for name in code_id]
ts_filenames = ["ts_" + name for name in code_id]
if opts.case == "desikan":
    ts_filenames = [filename + ".txt" for filename in ts_filenames]
else:
    ts_filenames = [filename + f"_{opts.case}.txt" for filename in ts_filenames]
	
# Read TS of these subjects and compute correlation matrix
list_corrs = []
for filename in tqdm.tqdm(ts_filenames):
    C_mat = np.corrcoef(np.loadtxt(opj(data_dir, f"{opts.case}" , filename)).T) # the original way, omiting the first point
    list_corrs.append(C_mat)
# Take upper triangular elements only
Y_cor = np.row_stack([squareform(mat, checks=False) for mat in list_corrs])
# Convert to Z-fisher
Y = np.arctanh(Y_cor)

filter_nans_subjects = ~np.any(np.isnan(Y), axis=1)

Y = Y[filter_nans_subjects]

subjects = subjects[filter_nans_subjects]
age = age[filter_nans_subjects]
sex = sex[filter_nans_subjects]
site = site[filter_nans_subjects]
group = group[filter_nans_subjects]
motion = motion[filter_nans_subjects]

# Compute dummies for sex and group
sex_dummies = pd.get_dummies(sex, drop_first=True).values
group_dummies = pd.get_dummies(group, drop_first=True).values

# X=Effects of Interest, C=Covariates
C = np.column_stack((sex_dummies, age, motion)) # Covariate effects
X = (group_dummies==0).astype(int).reshape(-1,1) # effects of interest, i.e. group varialbe

# Some constants
n_asd = (X==1).sum() # Number of ASD subjects
n_td = (X==0).sum()  # Number of TDC subjects
n_rois = C_mat.shape[1]  # Number of Rois
n_links = Y.shape[1] # Number of connectivity entries

print("ASD subjects: {}, TD subjects: {}, Number of rois: {}, Number of links: {}".format(n_asd, 
                                                                                          n_td, 
                                                                                          n_rois, n_links))

# %% HARMONISE
from pycombat import Combat
from sklearn.preprocessing import LabelEncoder

# Convert site labels to integers
le =  LabelEncoder()
b =le.fit_transform(site.values)

combat = Combat()
# Run combat, preserving  group (ASD/TDC) differences
Y_combat = combat.fit_transform(Y, b=b, X=X, C=None)

# %% DENOISE
from sklearn.linear_model import LinearRegression

del Y # This is to guarantee that I am not taking the original data

# Take only ASD subjects
Y_asd_combat = Y_combat[np.squeeze(X)==1,:]
print(Y_asd_combat.shape)

linReg = LinearRegression()
linReg.fit(C[np.squeeze(X)==1], Y_asd_combat)
Y_asd_combat_denoised = Y_asd_combat - C[np.squeeze(X)==1,:].dot(linReg.coef_.T)

# %% SAVE
CC_combat_asd_denoised = np.array([squareform(Y_asd_combat_denoised[ii, :]) + np.identity(n_rois) \
                                   for ii in range(n_asd)])
    
results_dir = f"/home/javi/Documentos/asd-subtyping-enrichment/results/{opts.case}"
Path(results_dir).mkdir(parents=True, exist_ok=True)

np.savez(opj(results_dir,
             'data_after_combat_motion_aggressive.npz'), 
         subjects=subjects,
         age = age, 
         sex=sex, 
         site=site, 
         group=group, 
         Y_combat = Y_combat,
         CC_combat = np.array([squareform(y) + np.eye(n_rois) for y in Y_combat]),
         motion = motion,
         CC_combat_asd_denoised = CC_combat_asd_denoised)

# Save connectivity matrices alone for clustering
savemat(opj(results_dir,'cors_combat_motion_aggressive.mat'), 
        {'CC_combat_asd_denoised': CC_combat_asd_denoised})
