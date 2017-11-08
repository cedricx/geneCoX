require('corrplot')
require('igraph')
require('R.matlab')
require('ggplot2')
setwd('~/Google Drive/GraduateSchool/Courses/BE566/Project1/')

load('./whole_brain.RData')
null<-readMat('./wb_null.mat')
brain_stats <- readMat('./wb_stats.mat')
coex_pwr <- abs(gene_coex_wb)^7
coex_graph<-graph_from_adjacency_matrix(coex_pwr,mode = c("undirected"), weighted = TRUE, diag = FALSE)
coex_comm<-cluster_fast_greedy(coex_graph, merges = TRUE, modularity = TRUE,membership = TRUE)
coex_mod <- modularity(coex_graph, membership =coex_comm$membership, weights = E(coex_graph)$weight)
coex_stg <- graph.strength(coex_graph,loops = FALSE)
coex_cent <- estimate_betweenness(coex_graph,cutoff=0,weights = E(coex_graph)$weight)

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

node_del <- order(coex_stg)[1:100]
coex_pwr_less <- coex_pwr[-node_del,-node_del]
coex_pwr_less[which(coex_pwr_less < quantile(coex_pwr_less,0.95))] <- 0
coex_graph_less <- graph_from_adjacency_matrix(coex_pwr_less,mode = c("undirected"), weighted = TRUE, diag = FALSE)
margin = .9
genes_final_na <- genes_final
genes_final_na[which(coex_stg < quantile(coex_stg,0.75))] <- ""
plot.igraph(coex_graph_less, vertex.label.cex = 0.5, vertex.label = genes_final_na[-node_del], 
            vertex.size = sqrt(coex_stg[-node_del]*5), vertex.label.color = "black",
            vertex.label.dist = 0, vertex.color = coex_comm$membership[-node_del],
            edge.width = E(coex_graph_less)$weight*20, edge.curved = TRUE,
            xlim=c(-margin,margin),ylim=c(-margin,margin))


#genes
gene_class <-read.csv('./pantherGeneList.csv',col.names = c("geneID","geneSymbol","geneName","Family","ProteinClass","Species"),header =FALSE)
genes_final_corrected <- genes_final
genes_final_corrected[202] <- "KLC1" #in Allan Brain data, KLC1 was mistakenly pointed to C14orf153
genes_final_class <- data.frame(gene_name = genes_final_corrected, gene_order = 1:331)
genes_final_class_merge <- merge(genes_final_class,gene_class, by.x = "gene_name", by.y = "geneSymbol")
genes_final_class_merge <- genes_final_class_merge[order(genes_final_class_merge$gene_order),]
genes_final_class_merge <- cbind(genes_final_class_merge,coex_comm$membership)

node_info <- data.frame(id = genes_final_class_merge$gene_order, label = genes_final_class_merge$gene_name, strength = coex_stg, centrality = brain_stats$wb.bt, 
                        com = com_data$CIU, chrom = genes_final_chrome$chrom)
write.csv(node_info,'./Figures/gephi_node.csv',quote = FALSE,row.names = FALSE)

#com_null
com_data <- readMat('./modularity_analysis.mat')
ave_null_modularity <- t(com_data$ave.null.mod)
ave_modularity <- unlist(com_data$wb.com[3,,])
pl.mod = data.frame(gamma = seq(from = 0.1, to = 10, by = 0.1), ave_mod = ave_null_modularity)
pl.mod$mod.lb = pl.mod$ave_mod - t(com_data$std.null.mod)*2
pl.mod$mod.up = pl.mod$ave_mod + t(com_data$std.null.mod)*2

#zoomed-out version
ggplot(pl.mod, aes(gamma)) + 
  geom_line(aes(y=ave_mod), colour="blue") + 
  geom_ribbon(aes(ymin=mod.lb, ymax=mod.up), alpha=0.2) +
  geom_line(aes(y=unlist(com_data$wb.com[3,,])), colour="red") +
  geom_line(aes(y=mod_diff), colour="orange") +
  ylab('Modularity')

#zoomed-in version
ggplot(pl.mod, aes(gamma)) + 
  geom_line(aes(y=ave_mod), colour="blue") + 
  geom_ribbon(aes(ymin=mod.lb, ymax=mod.up), alpha=0.2) +
  geom_line(aes(y=unlist(com_data$wb.com[3,,])), colour="red") +
  geom_line(aes(y=mod_diff), colour="orange") +
  coord_cartesian(xlim = c(0,1.5)) +
  ylab('Modularity')
mod_diff <- ave_modularity - ave_null_modularity 

#consensus partrition
consensus = com_data$consensus
consensus_ro = consensus
for (i in 1:dim(consensus)[1]) {
  com_i =  consensus[i,]
  com_sum = as.data.frame(sort(table(com_i)))
  for (j in 1:dim(com_sum)[1]){
    com_i[which(com_i == com_sum$com_i[j])] = 10+j
  }
  consensus_ro[i,] = com_i-10
}
  
consensus_ro = t(consensus_ro)

keycol=c('#EF5350','#4FC3F7','#AED581','#FFEB3B','#BA68C8',"#7986CB",'#FF6E40')#red,blue,green,yellow,purple, indego,orange
consensus_plot <-levelplot(consensus_ro, col.regions=rev(keycol),
                           at = seq(0,7,length.out = 8))

levelplot(com_data$D,par.settings = YlOrRdTheme())

ticks <- unname(table(com_data$CIU))
adj.keycol=c('212121','#FFF4B3','#FEE086','#FEBE59','#FD9740','#FC5C2E','#E51F1D','#BE0125','#800026')
levelplot(com_data$adj.pwr[order(com_data$CIU),order(com_data$CIU)], col.regions = adj.keycol,
          at = c(0,0.01,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8),xlab="",ylab = "", strip = F, contour = F, region= T,
          scales=list(x=list(tck = 0),
                      y=list(tck = 0)),
           panel = function(...){
            panel.levelplot(...)
            panel.abline(h = cumsum(ticks)+0.5,col="white")
            panel.abline(v = cumsum(ticks)+0.5,col="white")
          }
          )




