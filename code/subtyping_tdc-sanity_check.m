clear;clc;

working_dir = "/path/to/working/folder";

%%% Add functions to the path %%%%
addpath(genpath(strcat(working_dir, 'code')))

%%%% Load data %%%%%%
load(strcat(working_dir, "data/consensus_data_tdc_split_sanity_check_Jul22.mat"));

CC = CC_1_denoised; 
kind = 2:1:20;
nshuff= 1000;
gamma=1;
seed=0;


res_cluster_euclidean = cluster_consensus(CC, 'euclidean', kind, nshuff, gamma, seed);

save(strcat(working_dir, 'data/results_clustering_tdc_partition_k_0220_Jul22.mat'),...
    'res_cluster_euclidean')
