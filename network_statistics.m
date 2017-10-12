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

