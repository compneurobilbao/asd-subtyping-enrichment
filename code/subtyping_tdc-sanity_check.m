clear;clc;

working_dir = "/home/javi/Documentos/asd-subtyping-enrichment/";

%%% Add functions to the path %%%%
addpath(genpath(strcat(working_dir, 'code')))

%%%% Load data %%%%%%
load(strcat(working_dir, "results/desikan/consensus_data_tdc_split_sanity_check.mat"));

CC = CC_1_denoised; 
kind = 2:1:20;
nshuff= 1000;
gamma=1;
seed=0;


res_cluster_euclidean = cluster_consensus(CC, 'euclidean', kind, nshuff, gamma, seed);

save(strcat(working_dir, 'results/desikan/results_clustering_tdc_partition_k_0220_Jul22.mat'),...
    'res_cluster_euclidean')
