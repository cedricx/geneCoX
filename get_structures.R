working.dir <- '~/Google Drive/GraduateSchool/Courses/BE566'
setwd(working.dir)
source(paste0(working.dir,'/code/ABI_functions.R'))

ontology <- read.csv(paste0(working.dir,'/Project1/ABI/normalized_microarray_donor10021/Ontology.csv'))
ont.parse <- strsplit(as.character(ontology$structure_id_path),'/')

#Grey vs. White Matter vs. Spaces
struc.GM <- get_structures("4006") #grey matter
struc.WM <- get_structures("9218") #white matter
struc.SS <- get_structures("9352") #sulci & spaces

#Primitive brains
struc.Tel <- get_structures("4007") #telencephalon
struc.DiE <- get_structures("4391") #diencephalon
struc.MES <- get_structures("9001") #mesencephalon
struc.MET <- get_structures("4833") #metencephalon
struc.MY <- get_structures("9512") #myelencephalon

#Telecenphalon
struc.Cx <- get_structures("4008") #cerebral cortex (sub str. of telencephalon)
struc.CxN <- get_structures("4275") #cerebral nuclei (sub str. of telencephalon)

#Cerebral Cortex
struc.FL <- get_structures("4009") #frontal lobe (sub str. of Cx)
struc.Ins <- get_structures("4268") #insula (sub str. of Cx)
struc.LL <- get_structures("4219") #limbic lobe (sub str. of Cx)
struc.OL <- get_structures("4180") #occipital lobe (sub str. of Cx)
struc.PL <- get_structures("4084") #parietal lobe (sub str. of Cx)
struc.TL <- get_structures("4132") #temporal lobe (sub str. of Cx)

#for each donor construct gene coexpression network for each structure type
struc_list <- ls(pattern = 'struc.*',all.names = FALSE)
struc_list <- struc_list[-c(1,2)]
donor_list = c('9861','10021','12876','14380','15496','15697')
load(paste0(working.dir,'/Project1/gene_list_final.RData'))
genes_final <- c('C4A',genes_final)
genes_donor = NULL
# load processed data from each donor
for (donor in donor_list) {
  print(paste('Working on donor', donor))
  load(paste(working.dir,"/Project1/ABI/",donor,".RData",sep = ""))
  datExprIn <- datExpr[genes_final,]
  colnames(datExprIn) <- sampleInfo$structure_id
  genes_donor[[donor]] = list(donor = donor, exp = datExprIn, info = sampleInfo)
}

# calculate each structure for each donor
genes_donor_struc = NULL
for (donor in donor_list) {
  print(paste('Working on donor', donor))
  for (struc in struc_list) {
    print(paste('Working on struc', struc))
    struc_in_here <- colnames(genes_donor[[donor]]$exp) %in% as.character(get(struc)$id) 
    genes_donor_struc[[donor]][[struc]]$sample <- colnames(genes_donor[[donor]]$exp)[struc_in_here]
    if (length(which(struc_in_here == TRUE)) > 5){ #if this donor contain these samples, then calculate the coexpression network
    genes_donor_struc[[donor]][[struc]]$adj <- make_coex_graph(genes_final, genes_donor[[donor]]$exp[,struc_in_here],'pearson')
    } else { #if not, then skip
      genes_donor_struc[[donor]][[struc]]$adj <- NA
    }
  }
}

# calculate structures for all donor concatenated
genes_all = NULL
for (donor in donor_list){
  genes_all <- cbind(genes_all,genes_donor[[donor]]$exp)
}
genes_all_stru <- NULL
for (struc in struc_list) {
  print(paste('Working on struc', struc))
  struc_in_here <- colnames(genes_all) %in% as.character(get(struc)$id) 
  genes_all_stru[[struc]]$sample <- colnames(genes_all_stru)[struc_in_here]
  if (length(which(struc_in_here == TRUE)) > 5){ #if this donor contain these samples, then calculate the coexpression network
    genes_all_stru[[struc]]$adj <- make_coex_graph(genes_final, genes_all[,struc_in_here],'pearson')
  } else { #if not, then skip
    genes_all_stru[[struc]]$adj <- NA
  }
}
get_modularity<-function(coex_pwr){
coex_graph<-graph_from_adjacency_matrix(coex_pwr,mode = c("undirected"), weighted = TRUE, diag = FALSE)
coex_comm<-cluster_fast_greedy(coex_graph, merges = TRUE, modularity = TRUE,membership = TRUE)
coex_mod <- modularity(coex_graph, membership =coex_comm$membership, weights = E(coex_graph)$weight)
return(coex_mod)}

#ant-pos
sample_info_all <- NULL
genes_donor_AP <- NULL
for (donor in donor_list) {
  sample_info <- read.csv(paste0(working.dir,'/Project1/ABI/normalized_microarray_donor',donor,'/SampleAnnot.csv'))
  #sample_info_all <- rbind(sample_info_all,sample_info)
  genes_donor_AP[[donor]]<-get_AP_structures(sample_info)
}


genes_donor_AP_adj <- NULL
for (donor in donor_list[3:5]) {
  #struc_in_here <- lapply(genes_donor_AP[[donor]], function(AP_seg) colnames(genes_donor[[donor]]$exp) %in% as.character(AP_seg))
  #genes_donor_AP_exp <- lapply(struc_in_here, function(strc) genes_donor[[donor]]$exp[,strc])
  genes_donor_AP_exp <- lapply(genes_donor_AP[[donor]], function(AP_seg)  genes_donor[[donor]]$exp[,AP_seg])
  genes_donor_AP_adj[[donor]]$sample <- genes_donor_AP[[donor]]
  genes_donor_AP_adj[[donor]]$adj <- lapply(seq_along(genes_donor_AP_exp), function(i) {print(paste('working on AP segment #',i,'for donor',donor));
                                                                                    AP_seg_length <- length(genes_donor_AP_adj[[donor]]$sample[[i]]);
                                          if (AP_seg_length > 5) { out<-make_coex_graph(genes_final,genes_donor_AP_exp[[i]],'pearson') 
                                                                                  }  
                                                    else {
                                                    out <- NA;
                                                     print(paste('Sample',i, 'too short =',exp_dim,'for donor',donor))};
                                                                                    return(out)})
  
}
