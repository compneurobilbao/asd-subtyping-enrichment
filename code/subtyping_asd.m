clear;clc;

working_dir = "/home/javi/Documentos/asd-subtyping-enrichment/";

%%% Add functions to the path %%%%
addpath(genpath(strcat(working_dir, 'code')))

%%%% Load data %%%%%%
load(strcat(working_dir, "results/desikan/cors_combat_both_motion_aggressive.mat"));

% subset ASD
% asd_subset = group==1;
CC_combat_asd = CC_combat_asd_denoised; 
kind = 2:1:20;
nshuff= 1000;
gamma=1;
seed=0;


res_cluster_euclidean = cluster_consensus(CC_combat_asd, 'euclidean', kind, nshuff, gamma, seed);

save(strcat(working_dir, 'data/results_clustering_motion_aggressive_k_0220.mat'),...
    'res_cluster_euclidean')
