function res_cluster = cluster_consensus(CC, dist_type, kind, nshuff, gamma, seed)

    n_rois = size(CC, 2);
    n_subs = size(CC, 1);
    
    disp('starting cluster calculation with:')
    disp(strcat('n_rois=', num2str(n_rois)))
    disp(strcat('n_subs=', num2str(n_subs)))
    disp(strcat('dist_type=', dist_type))
    disp('kind =');disp(kind);
    disp(strcat('n_shuff=', num2str(nshuff)))
    disp(strcat('gamma=', num2str(gamma)))
    disp(strcat('seed=', num2str(seed)))
    
    D = zeros(n_subs, n_subs, n_rois);
    switch dist_type
        case 'pearson'
            for ii=1:n_rois
                D(:, :, ii) = sqrt(2*squareform(pdist(squeeze(CC(:, ii,1:n_rois~=ii)), 'correlation')));
            end 
        case 'spearman'
            for ii=1:n_rois
                D(:, :, ii) = sqrt(2*squareform(pdist(squeeze(CC(:, ii,1:n_rois~=ii)), dist_type)));
            end 
        otherwise 
            for ii=1:n_rois
                D(:, :, ii) = squareform(pdist(squeeze(CC(:, ii,1:n_rois~=ii)), dist_type));
            end 
    end

    tic
    [B, C, P] = consensus(D, kind, nshuff, gamma, seed);
    toc

    figure;
    imagesc(C)


    S_reps = {};
    Q_reps = [];

    rng(seed)
    n_reps = 100;
    for i=1:n_reps

        [S, Q] = genlouvain(B);

        S_reps{i} = S;
        Q_reps(i) = Q;
    end

    Rmat = zeros(n_reps);
    for ii=1:n_reps
        for jj=(ii+1):n_reps
        Rmat(ii,jj) = RandIndex(S_reps{ii}, S_reps{jj});
        Rmat(jj,ii) = Rmat(ii,jj);
        end
    end

    %we can see that the structure is always de same
    figure;
    imagesc(Rmat)

    Q_max_id = find(Q_reps==max(Q_reps));

    for i = 1:length(Q_max_id)

        [clus,~,ic] = unique(S_reps{Q_max_id(i)});
        a_counts = accumarray(ic,1);
        value_counts = [clus, a_counts];
        disp(value_counts)
    end

    clus_idxs = {};
    %clus_subj = {};
    for iclus = 1:length(unique(S_reps{Q_max_id(1)}))
        clus_idxs{iclus} = find(S_reps{Q_max_id(1)}==iclus);
        %clus_subj{iclus} = asd_subjects(clus_idxs{iclus});
    end

    S = S_reps{Q_max_id(1)};
    Q = Q_reps(Q_max_id(1));

    idxs= [];
    for ii=1:length(unique(S))
        idxs = [idxs;find(S==ii)];
    end

    figure;
    imagesc(C(idxs, idxs))
    title("ASD Consensus Matrix")
    xlabel("ASD subjects")
    ylabel("ASD subjects")
    colorbar

    res_cluster = struct;
    
    res_cluster.B = B;
    res_cluster.P = P;
    res_cluster.C = C;
    res_cluster.Q = Q;
    res_cluster.S = S;
    
end

