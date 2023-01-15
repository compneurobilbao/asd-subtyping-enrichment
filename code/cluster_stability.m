clear;clc;
working_dir = '/home/javi/Documentos/asd-subtyping-enrichment';

%%% Add functions to the path %%%%
addpath(genpath(strcat(working_dir, 'code')))

%%%% Load data %%%%%%
load(strcat(working_dir, '/results/desikan/',...
    'results_clustering_motion_aggressive_k_0220.mat'),...
    'res_cluster_euclidean');

B_obs = res_cluster_euclidean.B;
C_obs = res_cluster_euclidean.C;
S_obs = res_cluster_euclidean.S;
k_obs = full(sum(C_obs));
twom = sum(k_obs);
M_obs = res_cluster_euclidean.Q/twom;

n_asd = size(B_obs, 2);

nperms=1000;
rng(1234)
[~,idx]=bootstrp(nperms, [], 1:n_asd);

for ii=1:nperms
    B_perm = B_obs(idx(:,ii),idx(:,ii));
    [S_perm,Q_perm] = genlouvain(B_perm, 'verbose',0);
    
    C_perm = C_obs(idx(:,ii),idx(:,ii));
    k_perm = full(sum(C_perm));
    twom_perm = sum(k_perm);
    M_perm(ii) = Q_perm/twom_perm;
    S_perms(ii,:) = S_perm;
end
% Confidence intervals at alpha=0.05 in the Modularity:
alpha=0.05;
cis = quantile(M_perm, [alpha/2 1-alpha/2]);

disp(['observed modularity: ' num2str(M_obs)])
disp(['confidence intervals: [' num2str(cis(1)) ' , ' num2str(cis(2)) ']'])

% step 2: see stability of each cluster following XXX
n_clus = length(unique(S_obs));
stability_clus = NaN(n_clus, nperms);

for ii=1:nperms
    [C,IA,~] = unique(idx(:,ii));
    
    xprima = S_obs(C);
    delta = transpose(S_perms(ii, IA));
    n_clus_boot = length(unique(delta));
    
    for clus_id=1:n_clus
        if sum(xprima==clus_id)==0
            % Skip if this cluster was not presented in the bootstrapped
            % sample.
            continue
        end
        gamma = 0;
        for clus_boot_id=1:n_clus_boot
            temp = jaccard(xprima==clus_id, delta==clus_boot_id);
            if temp>gamma
                gamma = temp;
            end
        end
        stability_clus(clus_id, ii) = gamma;
    end
end

% Here number 2 in S_obs corresponded to our first reported subtype
% (Hypo), and number 1 to the second subtype (Hyper). This is because we
% later redefine these labels according to the number of subjects within
% each subgroup, with the first one being the largest.
disp(['Stability subtype 1 (Hypo) is: ', num2str(mean(stability_clus(2,:)))])
disp(['Stability subtype 2 (Hyper) is: ', num2str(mean(stability_clus(1,:)))])
