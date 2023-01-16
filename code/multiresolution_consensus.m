%% Using Consensus matrix
clear;clc;

working_dir = "/home/javi/Documentos/asd-subtyping-enrichment/";

%%% Add functions to the path %%%%
addpath(genpath(strcat(working_dir, 'code')))

load(strcat(working_dir, 'results/desikan/results_clustering_motion_aggressive_k_0220.mat'))

% Use the consensus clustering
A_consensus = res_cluster_euclidean.C;

S = eventSamples(A_consensus, 1000);
[Sc, Tree] = hierarchicalConsensus(S, 'SimilarityType', 'linkage');
C = coclassificationMatrix(S);

figure;
consensusPlot(C, Sc, Tree)

% Compute modularity for each partition
[Sall, thresholds] = allPartitions(Sc, Tree);

s_int=treeSort(C,  Sc, Tree);

save(strcat(working_dir, 'results/desikan/multi_resolution_consensus_motion_aggressive_k_0220.mat'),...
    'C','Sc','Tree', 'Sall', 'thresholds', 's_int')

f=figure;
f.Position = [100 100 1500 1300];
[ax_C, ax_H] =  consensusPlot(C, Sc, Tree);
ax_C.FontSize=20;
ax_H.FontSize=20;
xlabel(ax_C, "ASD Subjects", 'fontsize', 30)
yticklabels(ax_C, [])
%ylabel(ax_C, "ASD Subjects", 'fontsize', 30)
title(ax_C, 'Multiresolution Consensus Matrix', 'fontsize', 40)
title(ax_H, {['Multiresolution'] ['hierarchical'] ['Tree']}, 'fontsize', 20)
colorbar(ax_C, 'LOCATION', 'WestOutside')
set(gca, 'units', 'normalized'); %Just making sure it's normalized
Tight = get(gca, 'TightInset');
saveas(gcf, strcat(working_dir, 'plots/multires_consensusplot.png'))
saveas(gcf, strcat(working_dir, 'plots/multires_consensusplot.pdf'))
saveas(gcf, strcat(working_dir, 'plots/multires_consensusplot.svg'))
