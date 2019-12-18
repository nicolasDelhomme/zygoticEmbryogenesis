library(here)
library(RLinuxModules)
library(data.table)
module("load bioinfo-tools seidr-devel")
module("load bioinfo-tools InfoMap")
source(here("Rtoolbox/src/infomapTools.R"))
setwd("/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/")

projectDir <- "/mnt/picea/home/mstewart/Git/zygoticEmbryogenesis/data/seidr/"
#deDir <- paste0(projectDir, "de/")
networkDir <- paste0(projectDir, "network/")

#That is on the backbone-1-percent.sf file and the algorithm that we will use is the irp_score.
#Now we define those parameters, as well as the edge file and the names of the files that will come after Infomap.

#seidrExe <- ("/mnt/picea/home/bastian/Git/seidr-devel/build/seidr ")

seidrExe <- ("seidr ")
backboneFile <- paste0(projectDir, "results/backbone/backbone-1-percent.sf")
clusterFolder <- paste0(networkDir, "cluster/")
dir.create(clusterFolder, showWarnings = FALSE)
algo <- "irp_score"
edgeFile <- paste0(clusterFolder,"edgeList.txt")
edgeIndexFile <- paste0(clusterFolder,"edgeIndexList.txt")
treeFile <- paste0(clusterFolder,"edgeIndexList.tree")

# Infomap clustering
# First we use Seidr reheader function to drop not connected nodes and improve efficiency

system(paste(seidrExe, "reheader" , backboneFile), intern=TRUE)

#We execute seidr view to extract the information of source node, target node and the score of irp

headResult <- system(paste(seidrExe, "view", backboneFile, "-c -d $'\t' ", "| head -n 1"), intern=TRUE)
headResult <- unlist(strsplit(headResult, "\t"))
algoIndex <- grep("irp_score", headResult)
system(paste0(seidrExe, "view ", backboneFile, "  -d $'\t' | cut -f 1,2,",algoIndex, " > ",edgeFile), intern=TRUE)
system(paste0(seidrExe, "view ", backboneFile, " -N -d $'\t' | cut -f 1,2,",algoIndex," >",edgeIndexFile), intern=TRUE)

#We will execute Infomap with a markov-time value of 2, because that value matches our goal to 
#keeps between 50-60% of the genes in the top 20 clusters.

markovTime <- 0.556
system(paste("Infomap ", edgeIndexFile," -z --markov-time ", markovTime," ", clusterFolder))
infomapRes <- system(paste0(seidrExe, " resolve -s ", backboneFile, " ", treeFile), intern=TRUE)
infomapTable <-data.frame(do.call(rbind, strsplit(infomapRes, "\t")))
infomapTable <- prepareData(infomapTable)
infomapTable$Level1 <- infomapTable$P1
infomapTable$Level2 <- paste0(infomapTable$Level1,":",infomapTable$P2)
infomapTable$Level2 <- ifelse(infomapTable$Level2 %like% "NA", NA, infomapTable$Level2)  

#We can get some information about the quality of the clusters

print(paste("% of genes in the top 20 clusters in Level 1:", clusterQA(infomapTable)))
print(paste("% of genes in the top 20 clusters in Level 2:", clusterQA(infomapTable, level='Level2')))

#' MT 0.5, L1 96.94%, L2 48.92%, 
#' MT 0.55 L1 %96.94, L2 %56.77,
#' MT 0.555, L1 61.51%, L2 10.61% G20 20668, GT 33603
#' MT 0.556, L1 58.6%, L2 11.38% G20 19691, GT 33603 ##USING THIS ONE
#' MT 0.557, L1 58.22%, L2 13.00% G20 19562, GT 33603
#' MT 0.558, L1 63.24%, L2 11.87%
#' MT 0.559, L1 59.13%, L2 10.82%
#' MT 0.56, L1 58.54%, L2 10.17%
#' MT 0.561, L1 61.55%, L2 9.91%
#' MT 0.562, L1 59.85%, L2 11.6%
#' MT 0.57, L1 62.83%, L2 10.43%
#' MT 0.58, L1 63.4%, L2 11.48%
#' MT 0.59, L1 69.17%, L2 14.11%
#' MT 0.595, L1 66.25%, L2 13.12%
#' MT 0.6 L1 %66.63, L2 14.17
#' MT 0.61, L1 67.64%, L2 14.31%
#' MT 0.62, L1 73.14%, L2 46.2%
#' MT 0.64, L1 76.08%, L2 15.14%
#' MT 0.65, L1 %74.92, L2 %17.06



#We select those clusters that at least has 40 genes. Our goal is to obtain between 25-35 clusters

selectedClusters <- getClusterByMinSize(infomapTable, level = 'Level1',min=400)
topClusters <- getClusters(infomapTable, "Level1", numberOfClusters=selectedClusters)

topClusters
#Then we save the results of clustering for network visualization
here("results")
save4Cytoscape(topClusters, file=paste0(clusterFolder,"InfomapClusters.tsv"))
save4Cytoscape(topClusters, file=here("data/seidr/network/cluster", "InfomapClusters.tsv"))







clusterQA(infomapTable, level = "Level1")


topClusters[31]

  
source(here("UPSCb-common/src/R/gopher.R"))
enr_clusters <- lapply(topClusters[1:31], function(x){
  print(length(x))
  if(length(x) > 1)
    gopher(x, task = list('go', 'mapman','kegg'), background = NULL, url="pabies", alpha = 0.05)
  else
    NULL
})

enr_clusters
length(enr_clusters)

enr_clusters$Cluster3

#enriched the clusters that we have
enr2tsv(enr_clusters, file=paste0(clusterFolder,"enrichedClusters"))



#PLOTTING TREEMAPS OF Top Clusters
for(i in 1:length(enr_clusters)){
  
  x <- enr_clusters
  a <- length(x[[i]])
  dir <- "analysis/seidrclusters/"
  
  plotname <- names(x[i])
  
  if(a != 0){
    print(plotname)
    
    #plot and save go treemap
    if(is.null(nrow(x[[i]][[1]])) == FALSE){
      png(file=here(str_c(dir, "go_",plotname, ".png")),
          width=1000, height=700)
      plotEnrichedTreemap(x[[i]], enrichment = "go", namespace = "none")
      dev.off()
      print("go")
    }
    
    #plot and save mapman treemap
    if(is.null(nrow(x[[i]][[2]])) == FALSE){
      png(file=here(str_c(dir ,"mapman_",plotname, ".png")),
          width=1000, height=700)
      plotEnrichedTreemap(x[[i]], enrichment = "mapman", clusterColor = "#1B75BC")
      dev.off()
      print("map")
    }
    
  }
}






#extract first degree neighbours
source(here("Rtoolbox/src/getFDN.R"))
#edgelist, then genes
geneList <- read_csv(here("analysis/DE/FMG_Time_B2_vs_B(-1)genes.csv"))
geneList <- geneList$X1
geneList <- str_replace_all(geneList, "[.]1","")

edgeList <- read.table(here("data/seidr/network/cluster/edgeList.txt"),header=F,sep="\t")
edgeList <- edgeList[,1:2]

# data setup
# edgelist (2 columns, source and linked gene - 3rd column with value is removed)
# genes 


#subscript out of bounds error - why is that? Trying to access something that is out of the range of the matrix
#the error was that it was in a matrix
test <- getFDN(edgeList,geneList)
test

#cluster number 3 appears to be the yeastbody-like structure in cytoscape
#perform treemap of Cluster 3, which is interesting. Can do a treemap of all of them too, and archive them.
source(here("Rtoolbox/src/plotEnrichedTreemap.R"))

#Cluster 3, for GO and MAPMAN
plotEnrichedTreemap(enr_clusters$Cluster3, enrichment = "go", namespace = "none")
plotEnrichedTreemap(enr_clusters$Cluster3, enrichment = "mapman", clusterColor = "#1B75BC")

plotEnrichedTreemap(enr_clusters$Cluster30, enrichment = "go", namespace = "none")

enr_DE_ZF_F_B34 <- gopher(c(DE_FMG_B3_B4), task = list('go', 'mapman'), background = NULL, url="pabies", alpha = 0.05)
enr_DE_ZF_Z_B34 <- gopher(c(DE_ZE_B3_B4), task = list('go', 'mapman'), background = NULL, url="pabies", alpha = 0.05)

plotEnrichedTreemap(enr_DE_ZF_F_B34, enrichment = "go", namespace = "none")
plotEnrichedTreemap(enr_DE_ZF_F_B34, enrichment = "mapman", namespace = "none")
plotEnrichedTreemap(enr_DE_ZF_Z_B34, enrichment = "go", namespace = "none")
plotEnrichedTreemap(enr_DE_ZF_Z_B34, enrichment = "mapman", namespace = "none")


DE_FMG_B3_B4
DE_FMG_B7_B6
DE_ZE_B3_B4
DE_ZE_B7_B6












