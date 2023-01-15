#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 21 10:42:09 2022

@author: javi
"""

#%% PLOT MAPS
from nilearn.datasets import load_mni152_template, fetch_atlas_schaefer_2018
from nilearn.image import resample_to_img, load_img, index_img, new_img_like
from nilearn.plotting import plot_img

schaefer_img = fetch_atlas_schaefer_2018(n_rois=100)['maps']
desikan_img = index_img(load_img("/home/javi/Documentos/asd-subtyping-enrichment/data/Desikan_genetics_1mm.nii.gz"),

                        0)
cortex_ids = np.unique(schaefer_img.get_fdata())[1:]
subcortex_ids = np.array([35,36,37,38,39,40,41,76,77,78,79,80,81,82])

template_img = load_mni152_template(resolution=1)
schaefer_img = resample_to_img(schaefer_img, template_img, interpolation='nearest')
desikan_img = resample_to_img(desikan_img, template_img, interpolation='nearest')

schaefer_desikan_data = np.zeros(template_img.shape)

for roi_id in cortex_ids:
    schaefer_desikan_data[schaefer_img.get_fdata()==roi_id] = roi_id
    
for ii, roi_id in enumerate(subcortex_ids):
    schaefer_desikan_data[desikan_img.get_fdata()==roi_id] = 101 + ii
    
schaefer_desikan_img = new_img_like(template_img, schaefer_desikan_data)

res_Rsq_df_euclidean = pd.read_csv("/home/javi/Documentos/asd-subtyping-enrichment/results/schaefer_desikan/mdmr_rsquare_euclidean_schaefer_desikan.csv")
n_rois = 114
for ii in range(3):
     data = np.zeros_like(schaefer_desikan_img.get_fdata())
     for jj in range(n_rois):
         data[schaefer_desikan_img.get_fdata()==(jj+1)] = res_Rsq_df_euclidean.iloc[:, ii].values[jj]
    
     stat_img = new_img_like(atlas_img, data)
     plot_img(stat_img, 
              colorbar=True,
              threshold=1e-6,
              display_mode="mosaic",
              cmap=plt.cm.RdYlBu_r,
              bg_img = template_img,
              draw_cross=False, 
              black_bg=True,
              title=res_Rsq_df_euclidean.columns[ii])
    
     plt.savefig(opj("/home/javi/Documentos/asd-subtyping-enrichment/results/schaefer_desikan", 
                    f"R2_MDMR_map_{res_Rsq_df_euclidean.columns[ii]}.png"),
                 dpi=300)