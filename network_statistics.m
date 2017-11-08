[nRows,nCols] = size(adj);
adj(1:(nRows+1):nRows*nCols) = 0;
adj_ord = adj(order,order);
adj_pwr = abs(adj_ord).^7;

%null model
W0_wb = zeros([size(adj_pwr),1000]);
for i = 1:size(W0_wb,3)
    disp(i);
    W0_wb(:,:,i) = null_model_und_sign(adj_pwr);
end

%strength
[stg_wb_pos, stg_wb_neg] = strengths_und_sign(adj_ord);
hist(stg_wb_pos);hist(stg_wb_neg);

%betweeness centrality
wb_bt = betweenness_wei(1./adj_pwr);
null_bt = zeros(size(wb_bt,1),size(W0_wb,3));
for i = 1:size(W0_wb,3)
    disp(i);
    null_bt(:,i) = betweenness_wei(1./W0_wb(:,:,i));
end

%clustering coef
ccoef = clustering_coef_wu_sign(adj_pwr,3);
null_ccoef = zeros(size(ccoef,1),size(W0_wb,3));
for i = 1:size(W0_wb,3)
    disp(i);
    null_ccoef(:,i) = clustering_coef_wu_sign(W0_wb(:,:,i),3);
end

%global efficiency
wb_eglob = efficiency_wei(adj_pwr);
null_eglob = zeros(size(wb_eglob,1),size(W0_wb,3));
for i = 1:size(W0_wb,3)
    disp(i);
    null_eglob(:,i) = efficiency_wei(abs(W0_wb(:,:,i)));
end

%modularity
wb_com = [];
gammas = [0.1:0.1:10];
for i = 1:size(gammas,2)
    wb_com(i).gamma = gammas(i);
    [wb_com(i).community,wb_com(i).modularity] = community_louvain(adj_pwr,gammas(i));
end

[wb_com_g1,wb_mod_g1] = community_louvain(adj_pwr,1);
[wb_com_g2,wb_mod_g2] = community_louvain(adj_pwr,1);


wb_com_null = [];
for j = 1:size(W0_wb,3)
    disp(j)
    for i = 1:size(gammas,2)
        wb_com_null(j,i).gamma = gammas(i);
        [wb_com_null(j,i).community,wb_com_null(j,i).modularity] = community_louvain(W0_wb(:,:,i),gammas(i));
    end
end


%null model
W0_wb = zeros([size(adj_pwr),1000]);
for i = 1:size(W0_wb,3)
    disp(i);
    W0_wb(:,:,i) = abs(null_model_und_sign(adj_ord,1,0.5)).^7;
end

%modularity in the null
wb_com_null = [];
for j = 1:size(W0_wb,3) %1000 permutations
    disp(j)
    for i = 1:size(gammas,2) %100 gammas
        wb_com_null(j,i).gamma = gammas(i);
        [wb_com_null(j,i).community,wb_com_null(j,i).modularity] = community_louvain(W0_wb(:,:,j),gammas(i));
    end
end

null_modularity = reshape([wb_com_null(:).modularity;],[1000,100]);
ave_null_mod = mean(null_modularity,1);
std_null_mod = std(null_modularity,1);


%community consensus analysis
wb_com_gamma = [];
for i = 1:100
    [wb_com_gamma(i).com, wb_com_gamma(i).mod] = community_louvain(adj_pwr,0.9);
end
consensus = [wb_com_gamma(:).com;];


D = agreement_weighted(consensus,[wb_com_gamma.mod]);
CIU = consensus_und(D,0.1,100);

save('modularity_analysis.mat','ave_null_mod','std_null_mod','consensus','wb_com','D','CIU','adj_pwr')




