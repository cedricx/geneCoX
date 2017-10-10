library(WGCNA)
library(RJSONIO)
setwd('~/Google Drive/GraduateSchool/Courses/BE566/Project1/')
source('./WGCNA/ABI_functions.R')

##################################
######## PreProcess Data #########
##################################
donor_list = c('9861','10021','12876','14380','15496','15697')
for (donor in donor_list){
  inputFile = paste('./ABI/normalized_microarray_donor',donor,sep="")
  outputFile = paste('./ABI/',donor,'.RData',sep="")
  preprocessData(inputFile,outputFile)
}

##################################
####### Compile Gene List ########
##################################
nature108genes <- read.csv ('~/Google Drive/GraduateSchool/Courses/BE566/Project1/nature108genes.csv')
genes0 <- unlist(nature108genes[,6:31])
genes <- as.character(unname(genes0[-which(genes0 =="")]))

gene_incl_idx <- which(genes %in% probeInfo$gene_symbol)
gene_in <- genes[gene_incl_idx]
gene_out <- genes[-gene_incl_idx]
write.csv(gene_out,file = './gene_out.csv')

gene_out_dict <- read.csv('./gene_out_dictionary.csv')
genes_trans <- genes #translated
genes_trans[genes %in% as.character(gene_out_dict$nature)] <- as.character(gene_out_dict$ALAN)
genes335 <- unname(genes_trans[-which(genes_trans =="")])

save(genes335,file = './gene_list_335.RData')

##############################################
##### Load and Compile PreprocessedData ######
##############################################
load('./gene_list_335.RData')
genes_donor = NULL
for (donor in donor_list) {
  load(paste("./ABI/",donor,".RData",sep = ""))
  datExprIn <- datExpr[genes335,]
  colnames(datExprIn) <- sampleInfo$structure_id
  genes_donor[[donor]] = list(donor = donor, exp = datExprIn, info = sampleInfo)
}
save(genes_donor, file = './genes_donor.RData')

##############################################
##### Construct Correlation Matrix ######
##############################################
load('./genes_donor.RData')

#whole brain for all donors
genes_all = NULL
for (donor in donor_list){
  genes_all <- cbind(genes_all,genes_donor[[donor]]$exp)
}

gene_coex_wb <- make_coex_graph(genes335,genes_all,"pearson")
  
wb_heat <- heatmap(gene_coex_wb)
p.gene_coex_wb<-levelplot(gene_coex_wb[wb_heat$rowInd,wb_heat$rowInd],par.settings = BuRdTheme(),
          xlab = "", ylab = "", main = "SCZ 335",
          scales=list(x=list(at=1:length(genes335),labels=genes335[wb_heat$rowInd],rot=90, tck = 0, cex =0.2),
                      y=list(at=1:length(genes335),labels=genes335[wb_heat$rowInd], tck = 0,cex =0.2)))

save(genes335,gene_coex_wb,p.gene_coex_wb,wb_heat,file = './whole_brain.RData')
genes335.ro <-genes335[wb_heat$rowInd]
writeMat('./whole_brain.mat', genes = genes335, adj = gene_coex_wb)


#for each donor
gene_coex_single = NULL
for (donor in donor_list){
  mat <- make_coex_graph(genes335,genes_donor[[donor]]$exp,"pearson")
  plot <- plot_coex_graph(donor,genes335,mat,wb_heat$rowInd)
  gene_coex_single[[donor]] = list(mat = mat, plot = plot)
}
save(gene_coex_single,wb_heat,file = './single_brains.RData')

#for each donor
structure_IDs = NULL
structure_names = NULL
for (donor in donor_list){
  structure_IDs <- c(structure_IDs,genes_donor[[donor]]$info$structure_id)
  structure_names <- c(structure_names,genes_donor[[donor]]$info$structure_name)
}
structures <- NULL
structures <- data.frame(ids = unique(structure_IDs), names = unique(structure_names))
structures <- structures[order(structures$ids),]
write.csv(structures, file = './structures_list.csv',row.names = FALSE)



powers = 1:20
soft_in = pickSoftThreshold(t(datExprIn),networkType="signed",powerVector=powers)
cex1 = 0.9
plot(soft_in$fitIndices[,1], -sign(soft_in$fitIndices[,3])*soft_in$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(soft_in$fitIndices[,1], -sign(soft_in$fitIndices[,3])*soft_in$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.80,col="red")
thresholdPower_in = 15

net_in = blockwiseModules(t(datExprIn),power=thresholdPower_in,networkType="signed")
dissTOM = 1-TOMsimilarityFromExpr(datExprIn, power = thresholdPower_in)
plotTOM = dissTOM^thresholdPower_in
diag(plotTOM) = NA


adj = adjacency(t(datExprIn),type = "signed", power = thresholdPower_in)
diag(adj) = 0
adj_order = adj[order(net_in$colors),order(net_in$colors)]
levelplot(adj_order[1:60,1:60], par.settings=BuRdTheme())
