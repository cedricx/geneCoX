require('corrplot')
require('igraph')
load('./whole_brain.RData')
null<-readMat('./wb_null.mat')
coex_pwr <- abs(gene_coex_wb)^7
coex_graph<-graph_from_adjacency_matrix(coex_pwr,mode = c("undirected"), weighted = TRUE, diag = FALSE)
coex_comm<-cluster_fast_greedy(coex_graph, merges = TRUE, modularity = TRUE,membership = TRUE)
coex_mod <- modularity(coex_graph, membership =coex_comm$membership, weights = E(coex_graph)$weight)
coex_stg <- graph.strength(coex_graph,loops = FALSE)

#visualization
corrplot(coex_pwr[order(coex_comm$membership),order(coex_comm$membership)],method = 'shade',tl.pos = 'n')
levelplot(coex_pwr[order(coex_comm$membership),order(coex_comm$membership)])

#null model
null_coex_graph = NULL;
null_comm = NULL;
null_coex_mod = NULL;
for (i in 1:dim(null$W0.wb)[3]){
  print(paste("Processing",i,"th null graph", sep=""))
  null_coex_graph[[i]] <- graph_from_adjacency_matrix(null$W0.wb[,,i],mode = c("undirected"), weighted = TRUE, diag = FALSE)
  null_comm[[i]] <- cluster_fast_greedy(null_coex_graph[[i]], merges = TRUE, modularity = TRUE,membership = TRUE)
  null_coex_mod[[i]] <- modularity(null_coex_graph[[i]], membership =null_comm[[i]]$membership, weights = E(null_coex_graph[[i]])$weight)
}

node_del <- order(coex_stg)[1:50]
coex_graph_less <- graph_from_adjacency_matrix(coex_pwr[-node_del,-node_del],mode = c("undirected"), weighted = TRUE, diag = FALSE)
margin = .8
genes_final_na <- genes_final
genes_final_na[which(coex_stg < quantile(coex_stg,0.75))] <- ""
plot.igraph(coex_graph_less, vertex.label.cex = 0.5, vertex.label = genes_final_na[-node_del], 
            vertex.size = sqrt(coex_stg[-node_del]*5), vertex.label.color = "black",
            vertex.label.dist = 0, vertex.color = coex_comm$membership[-node_del],
            edge.width = E(coex_graph_less)$weight*20, edge.curved = TRUE,
            xlim=c(-margin,margin),ylim=c(-margin,margin))
