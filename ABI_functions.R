preprocessData <- function(donorDir, outputFileName) {
  datExpr <- read.csv(sprintf("%s/MicroarrayExpression.csv",donorDir),header=FALSE,row.names=1)
  probeInfo <- read.csv(sprintf("%s/Probes.csv",donorDir),header=TRUE)
  sampleInfo <- read.csv(sprintf("%s/SampleAnnot.csv",donorDir),header=TRUE)
  
  validProbes = which(probeInfo$entrez_id>=0,arr.ind=TRUE)
  datExpr <- datExpr[validProbes,]
  probeInfo <- probeInfo[validProbes,]
  
  probeNames = probeInfo[which(probeInfo == rownames(datExpr)),"gene_symbol"]
  probeIDs = rownames(datExpr)
  datCollapsed = collapseRows(datExpr, probeNames, probeIDs)
  
  datExpr <- datCollapsed$datETcollapsed
  probeInfo <- probeInfo[which(datCollapsed$selectedRow == TRUE,arr.ind=TRUE),]
  
  sampleInfo <- downloadSampleProperties(sampleInfo)
  save(datExpr,probeInfo,sampleInfo,file=outputFileName)
}

getApiPath <- function() {
  return("http://api.brain-map.org")
}

apiQuery <- function(query) {
  totalRows <- -1
  rowsPerPage <- 2000
  startRow <- 0
  queryFormat <- "%s/api/v2/data/%s&startRow=%d&numRows=%d"
  
  while (totalRows < 0 || startRow < totalRows) {
    queryString <- sprintf(queryFormat, getApiPath(), query, startRow, rowsPerPage)
    print(queryString)
    
    resultString <- readLines(queryString)
    resultJSON <- fromJSON(resultString)
    
    if (totalRows < 0) {
      totalRows <- as.integer(resultJSON$total_rows)
      output <- resultJSON$msg
    } else {
      output <- c(output,resultJSON$msg)
    }
    
    startRow <- startRow + rowsPerPage
  }
  
  return(output)
}

downloadStructure <- function(structureID) {
  queryString <- sprintf("%s/api/v2/data/Structure/%s.json", getApiPath(), structureID)
  print(queryString)
  resultString <- readLines(queryString,warn = FALSE)
  resultJSON <- fromJSON(resultString)
  
  return(resultJSON$msg)
}

downloadSampleProperties <- function(sampleInfo) {
  structureIDs = sampleInfo$structure_id
  numSamples <- length(structureIDs)
  
  structureColors <- array("#000000FF",numSamples)
  structureOrders <- array(0,numSamples)
  
  for (i in 1:numSamples) {
    s = downloadStructure(structureIDs[i])
    structureOrders[i] <- s[[1]]$graph_order
    structureColors[i] <- paste("#",s$color_hex_triplet,sep="")
  }
  
  sampleInfo[,"order"] <- structureOrders
  sampleInfo[,"color"] <- structureColors
  return(sampleInfo)
}

make_coex_graph <- function(geneList, coex_data, method){
  coex_graph <- matrix(NA,length(geneList),length(geneList))
  for (gene_i in geneList) {
    for (gene_j in geneList) {
      gene_cor <- cor(coex_data[gene_i,],coex_data[gene_j,],method = method)
      coex_graph[which(geneList == gene_i),which(geneList == gene_j)] <- gene_cor
    }
  }
  return(coex_graph)
}

plot_coex_graph <- function(donor,geneList,coex_graph,idx){
  plot<-levelplot(coex_graph[idx,idx],par.settings = BuRdTheme(),
                            xlab = "", ylab = "", main = donor,
                            scales=list(x=list(at=1:length(geneList),labels=geneList[idx],rot=90, tck = 0, cex =0.2),
                                        y=list(at=1:length(geneList),labels=geneList[idx], tck = 0,cex =0.2)))
  return(plot)
}


