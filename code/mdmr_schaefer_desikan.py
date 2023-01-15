#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 11:20:11 2022

@author: javi
"""
# %% IMPORTS
import os
import numpy as np
import pandas as pd
import argparse
from scipy.io import loadmat
import json
import matplotlib.pylab as plt
import numpy as np
from joblib import Parallel, delayed
from scipy.spatial.distance import pdist, squareform
from os.path import join as opj
from nilearn import image, plotting
from scipy.stats import ttest_ind

from mdmr import MDMR

# %% PARSE CASE
parser = argparse.ArgumentParser(description='Harmonise and denoise ASD before clustering')
parser.add_argument('--case',
                type=str,
                choices = ["desikan", "schaefer", "schaefer_desikan"],
                help="Name for the output directory")

opts = parser.parse_args()

# %% READ DATA
data_dir = "/home/javi/Documentos/asd-subtyping-enrichment/data"
results_dir = f"/home/javi/Documentos/asd-subtyping-enrichment/results/{opts.case}"

data_after_combat = np.load(opj(results_dir, 
                                'data_after_combat_motion_aggressive.npz'))

# harmonised FC data
Y_combat = data_after_combat['Y_combat']
print(Y_combat.shape)
n_rois = squareform(Y_combat[0,:], checks=False).shape[1]
n_subjects = Y_combat.shape[0]

# Convert these to connectivity matrices
CC_adjust = np.array(
    [squareform(Y_combat[ii, :], checks=False) + np.eye(n_rois) 
     for ii in range(n_subjects)]
    )
print(CC_adjust.shape)

# Covariates:
age = data_after_combat['age']
sex = data_after_combat['sex']
group = data_after_combat['group']
motion = data_after_combat['motion']

# Get demographics data to extract FIQ to be used as another covariate
subjects = data_after_combat['subjects']

asd_subjects = subjects[group==1]
# Load demo data
demo = pd.read_excel(opj(data_dir, "demo_final_FD03.xlsx"))
# Rearrange demographic data according to subjects in the data
demo = pd.merge(pd.DataFrame({'SUB_ID': subjects}), demo)

# %% ADD CLUSTERING INFO
# Add subtype IDs to demo data
def return_subject_clusters(S, threshold = 20):
    # Function to label the subtypes. Those below 20 were considered outliers
    cluster_id, cluster_counts = np.unique(S, return_counts=True)
    cluster_order = cluster_id[np.argsort(cluster_counts)[::-1]]
    
    clus_idxs = [np.where(S==ii)[0] for ii in cluster_order]
    
    clus_subj = [asd_subjects[idxs] for idxs in clus_idxs]
    
    cluster_labels = [[ii+1]*len(subjects) if len(subjects) > threshold else ['outlier']*len(subjects) \
                     for ii, subjects in enumerate(clus_subj)]
    
    return {'clus_subj': clus_subj, 'cluster_labels': cluster_labels, 'clus_idxs': clus_idxs}

case = 'res_cluster_euclidean'
cons_res=loadmat(opj(results_dir, 'results_clustering_motion_aggressive_k_0220.mat'), 
                 struct_as_record=False, 
                 squeeze_me=True)[case]

S_euclidean = cons_res.S
C_euclidean = cons_res.C
subjects_clusters_euclidean = return_subject_clusters(S_euclidean)
# Modularity from GenLouvain Matlab implementation
print(f"Modularity = {cons_res.Q/np.sum(cons_res.C)}")

demo['clus_id_euclidean'] = pd.NA

demo.loc[group==2, 'clus_id_euclidean'] = 0

for sub_clus, clus_label in zip(subjects_clusters_euclidean['clus_subj'], 
                                subjects_clusters_euclidean['cluster_labels']):
    demo.loc[demo.SUB_ID.isin(sub_clus), 'clus_id_euclidean'] = clus_label
    
#demo_data.clus_id_euclidean.value_counts()

avg_conn_pos = np.array([Y_combat[ii,:][Y_combat[ii,:]>0].mean() \
                for ii in range(Y_combat.shape[0])])
    
print("Group-average connectivity, subtype 1", 
      np.mean(avg_conn_pos[(demo.clus_id_euclidean==1).to_numpy()]))

print("Group-Average connectivity, subtype 2", 
      np.mean(avg_conn_pos[(demo.clus_id_euclidean==2).to_numpy()]))
print(" Testing these differences:", 
      ttest_ind(avg_conn_pos[(demo.clus_id_euclidean==1).to_numpy()], 
                avg_conn_pos[(demo.clus_id_euclidean==2).to_numpy()])
      )

# %% PREPARE DATA FOR MDMR
# Covariates: Age, sex, FIQ, FD and DVARS
covariates = ['AGE_AT_SCAN', 'SEX',  'FIQ', "FD_mean", "DVARS_mean"]
# input variables for MDMR: group and subtyping id (variables of interest) + covariates
variables = ['clus_id_euclidean', "DX_GROUP"] + covariates
# Discard any subject with missing info in some of these variables (normally only missing info in FIQ)
to_take = ~np.any(demo.loc[:, variables].isna(), axis=1).values

# Construct input matrix as data frame
XX_df = demo.loc[to_take, variables]
# Square (harmonized) connectivity matrices
YY = CC_adjust[to_take, :, :]

# Make some conversions for categorical data

XX_df['clus_id_euclidean'] = XX_df['clus_id_euclidean'].astype(str)
XX_df['SEX'] = XX_df['SEX'].astype(str)
XX_df['DX_GROUP'] = (XX_df['DX_GROUP']==1).astype(int).astype(str)

print("Subtypes (0==TDC)")
print(XX_df['clus_id_euclidean'].value_counts())
print("")

print("Groups")
print(XX_df['DX_GROUP'].value_counts())

X_df = XX_df.loc[:, ["clus_id_euclidean"] + covariates].copy()

# %% RUN MDMR

# Empty data frames for results
res_F_df_euclidean = pd.DataFrame({})
res_Rsq_df_euclidean = pd.DataFrame({})

# Take only subtype 1 and 2, which are the major ones
for jj in [1, 2]:
    print("doing cluster %d" % jj)
    list_Fs = []
    list_Rs = []
    # 0 is TDC, we compare with one of the subtypes
    cond = (X_df.loc[:, "clus_id_euclidean"] =='0') | (X_df.loc[:, "clus_id_euclidean"] ==str(jj))
    
    XX = X_df.loc[cond,:]
    for ii in range(YY.shape[1]):
        
        # Do not consider autocorrelation entries, which are 1
        mask = np.arange(n_rois)!=ii
        
        # Construct matrix distance
        D = squareform(pdist(YY[cond, ii, :][:, mask], 'euclidean'))
        mdmr = MDMR(verbose=0)
        mdmr.fit(XX, D)
        list_Fs.append(mdmr.F_[0])
        list_Rs.append(mdmr.r2_[0])
        
    res_F_df_euclidean['clus_' + str(jj)] = list_Fs
    res_Rsq_df_euclidean['clus_' + str(jj)] = list_Rs
    
# mdmr.r2_[0] above corresponds to our variable of interest, i.e. subtype ID versus TDC. Let's see it
print(mdmr.summary(), mdmr.r2_[0])

# Adding whole-group case
X_df =  XX_df.loc[:, ["DX_GROUP"] + covariates].copy()
print("")
print("doing whole group-case")
list_Fs = []
list_Rs = []

XX = X_df.copy()
for ii in range(YY.shape[1]):
    
    # Do not consider autocorrelation entries, which are 1
    mask = np.arange(n_rois)!=ii
    # Construct matrix distance
    D = squareform(pdist(YY[:, ii, :][:, mask], 'euclidean'))
    mdmr = MDMR(verbose=0)
    mdmr.fit(XX, D)
    list_Fs.append(mdmr.F_[0])
    list_Rs.append(mdmr.r2_[0])

res_F_df_euclidean.append(list_Fs)
res_Rsq_df_euclidean['whole'] = list_Rs

#%% SAVE RESULTS
res_F_df_euclidean.to_csv(opj(results_dir,
                              f"mdmr_stats_euclidean_{opts.case}.csv"), 
                          index=False)
res_Rsq_df_euclidean.to_csv(opj(results_dir, 
                                f"mdmr_rsquare_euclidean_{opts.case}.csv"), 
                            index=False)

# #%% PLOT MAPS
# from nilearn.datasets import load_mni152_template, fetch_atlas_schaefer_2018
# from nilearn.image import resample_to_img

# atlas_img = fetch_atlas_schaefer_2018(n_rois=100)['maps']

# template_img = load_mni152_template(resolution=1)
# atlas_img = resample_to_img(atlas_img, template_img, interpolation='nearest')

# for ii in range(3):
#     data = np.zeros_like(atlas_img.get_fdata())
#     for jj in range(n_rois):
#         data[atlas_img.get_fdata()==(jj+1)] = res_Rsq_df_euclidean.iloc[:, ii].values[jj]
    
#     stat_img = image.new_img_like(atlas_img, data)
#     plotting.plot_img(stat_img, 
#                       colorbar=True,
#                       threshold=1e-6,
#                       display_mode="mosaic",
#                       cmap=plt.cm.RdYlBu_r,
#                       bg_img = template_img,
#                       draw_cross=False, 
#                       black_bg=True,
#                       title=res_Rsq_df_euclidean.columns[ii])
    
#     plt.savefig(opj(results_dir, 
#                    f"R2_MDMR_map_{res_Rsq_df_euclidean.columns[ii]}.png"),
#                 dpi=300)
    
#%% SANITY CHECK USING R
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
ro.numpy2ri.activate() 
R = ro.r
R('library(MDMR)')

var_interest = "clus_id_euclidean"
res_Rsq_df_euclidean_rpy2 = pd.DataFrame({})

X_df = XX_df.loc[:, [var_interest] + covariates].copy()
n_clusters = len(pd.to_numeric(X_df.loc[:, var_interest], errors='coerce').dropna().unique())
print(pd.to_numeric(X_df.loc[:, var_interest], errors='coerce').dropna().value_counts())

for jj in [1, 2]:
    print("doing cluster %d" % jj)
    list_Rs_rpy2 = []
    cond = (X_df.clus_id_euclidean =='0') | (X_df.clus_id_euclidean ==str(jj))
    
    XX = X_df.loc[cond,:]
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_from_pd_df = ro.conversion.py2rpy(XX)
    R.assign('r_from_pd_df', r_from_pd_df)
    R('Xm <- model.matrix(~., data=r_from_pd_df)')
    
    for ii in range(YY.shape[1]):
            mask = np.arange(n_rois)!=ii
            D = squareform(pdist(YY[cond, ii, :][:, mask], 'euclidean'))
            R.assign('D', D)    
            
            list_Rs_rpy2.append(np.asarray(R('as.vector(mdmr(Xm, D, perm.p=FALSE)$pr.sq)'))[1][0])
    res_Rsq_df_euclidean_rpy2['clus_' + str(jj)] = list_Rs_rpy2
    
# Adding whole-group case
X_df =  XX_df.loc[:, ["DX_GROUP"] + covariates].copy()
XX = X_df.copy()
print("doing whole group-case")
list_Rs_rpy2 = []
with localconverter(ro.default_converter + pandas2ri.converter):
    r_from_pd_df = ro.conversion.py2rpy(XX)
R.assign('r_from_pd_df', r_from_pd_df)
R('Xm <- model.matrix(~., data=r_from_pd_df)')

for ii in range(YY.shape[1]):
        mask = np.arange(n_rois)!=ii
        D = squareform(pdist(YY[:, ii, :][:, mask], 'euclidean'))
        R.assign('D', D)    

        list_Rs_rpy2.append(np.asarray(R('as.vector(mdmr(Xm, D, perm.p=FALSE)$pr.sq)'))[1][0])
res_Rsq_df_euclidean_rpy2['whole'] = list_Rs_rpy2
    
print("For subtype 1, MDMR maps using our implementation and R are similar? ", 
      np.allclose(res_Rsq_df_euclidean.clus_1.to_numpy(), res_Rsq_df_euclidean_rpy2.clus_1.to_numpy()))
print("For subtype 2, MDMR maps using our implementation and R are similar?", 
      np.allclose(res_Rsq_df_euclidean.clus_2.to_numpy(), res_Rsq_df_euclidean_rpy2.clus_2.to_numpy()))
print("For the full ASD, MDMR maps using our implementation and R are similar? ", 
      np.allclose(res_Rsq_df_euclidean.whole.to_numpy(), res_Rsq_df_euclidean_rpy2.whole.to_numpy()))

