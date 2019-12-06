#Identification of TFs that may appear in Differential Expression

#Get List of TFs
TF_Baselist <- read.table(here("doc/TF_gene_names_and_families_names_picea_1841_HC_MC_LC"),header=F,sep="\t")
TF_Baselist <- as.character(TF_Baselist[,1])

testgenes1 <- genes_PC1_50_names
testgenes2 <- genes_PC2_50_names
testgenes3 <- genes_PC3_50_names

setdiff(testgenes1,testgenes2)
setdiff(testgenes1,testgenes3)
setdiff(testgenes2,testgenes3)

testgenes1 <- str_replace(testgenes1,"[.]1","")
testgenes2 <- str_replace(testgenes2,"[.]1","")
testgenes3 <- str_replace(testgenes3,"[.]1","")


tg1_res <- intersect(testgenes1,TF_Baselist)
tg2_res <- intersect(testgenes2,TF_Baselist)
tg3_res <- intersect(testgenes3,TF_Baselist)

tg1_res #81
tg2_res #105
tg3_res #65



setdiff(tg1_res,union(tg3_res,tg2_res)) #44
setdiff(tg2_res,union(tg1_res,tg3_res)) #62
setdiff(tg3_res,union(tg1_res,tg2_res)) #30


write.csv(tg1_res, file=here("analysis/ZE_FMG_PC1_TFs.tsv"))
write.csv(tg2_res, file=here("analysis/ZE_FMG_PC2_TFs.tsv"))
write.csv(tg3_res, file=here("analysis/ZE_FMG_PC3_TFs.tsv"))




#loop through and gather all the DE files into a single variable


FL_G_TE <- list.files(here("analysis/DE/ZF_TissueEffect"), 
                          recursive = TRUE, 
                          pattern = "",
                          full.names = TRUE)

FL_G_TE <- str_subset(FL_G_TE, "genes.csv")
FL_G_TE <- str_subset(FL_G_TE, "up", negate = TRUE)
FL_G_TE <- str_subset(FL_G_TE, "down", negate = TRUE)

FL_G_TE <- str_sort(FL_G_TE, numeric = TRUE)

de_genes_fmg_vs_ze <- sapply(1:length(FL_G_TE),function(i){
  
  GeneL <- read.csv(FL_G_TE[i])$X
  GeneL <- as.character(GeneL)
  GeneL <- str_replace(GeneL,"[.]1","")
  
  print(GeneL)
  return(GeneL)
  
})

for(i in 1:10){
  print(str_c("B",i))
  print(intersect(de_genes_fmg_vs_ze[i],TF_Baselist))
}





#import FMG vs Time genes
{
FL_G_TimeF <- list.files(here("analysis/DE/ZF_FMG_Time"), 
                      recursive = TRUE, 
                      pattern = "",
                      full.names = TRUE)

FL_G_TimeF <- str_subset(FL_G_TimeF, "genes.csv")
FL_G_TimeF <- str_subset(FL_G_TimeF, "up", negate = TRUE)
FL_G_TimeF <- str_subset(FL_G_TimeF, "down", negate = TRUE)

FL_G_TimeF <- str_sort(FL_G_TimeF, numeric = TRUE)
}
#clip off the .1 off the gene ids
de_genes_zf_fmg_vs_time <- sapply(1:length(FL_G_TimeF),function(i){
  
  GeneL <- read.csv(FL_G_TimeF[i])$X
  GeneL <- as.character(GeneL)
  GeneL <- str_replace(GeneL,"[.]1","")
  
  print(GeneL)
  return(GeneL)
  
})
#compare them to TF baselist
de_genes_zf_fmg_vs_time_TFs <- sapply(1:length(FL_G_TimeF),function(i){
  print(str_c("B",i))
  x <- intersect(unlist(de_genes_zf_fmg_vs_time[i]),TF_Baselist)
  print(x)
  return(x)
})

#import ZE vs Time genes
{
  FL_G_TimeZ <- list.files(here("analysis/DE/ZF_ZE_Time"), 
                           recursive = TRUE, 
                           pattern = "",
                           full.names = TRUE)
  
  FL_G_TimeZ <- str_subset(FL_G_TimeZ, "genes.csv")
  FL_G_TimeZ <- str_subset(FL_G_TimeZ, "up", negate = TRUE)
  FL_G_TimeZ <- str_subset(FL_G_TimeZ, "down", negate = TRUE)
  
  FL_G_TimeZ <- str_sort(FL_G_TimeZ, numeric = TRUE)
}
#clip off the .1 off the gene ids
de_genes_zf_ze_vs_time <- sapply(1:length(FL_G_TimeZ),function(i){
  
  GeneL <- read.csv(FL_G_TimeZ[i])$X
  GeneL <- as.character(GeneL)
  GeneL <- str_replace(GeneL,"[.]1","")
  
  print(GeneL)
  return(GeneL)
  
})
#compare them to TF baselist
de_genes_zf_ze_vs_time_TFs <- sapply(1:length(FL_G_TimeZ),function(i){
  print(str_c("B",i))
  x <- intersect(unlist(de_genes_zf_ze_vs_time[i]),TF_Baselist)
  print(x)
  return(x)
})

#import ZE vs FMG genes
{
  FL_G_Tissue <- list.files(here("analysis/DE/ZF_TissueEffect"), 
                           recursive = TRUE, 
                           pattern = "",
                           full.names = TRUE)
  
  FL_G_Tissue <- str_subset(FL_G_Tissue, "genes.csv")
  FL_G_Tissue <- str_subset(FL_G_Tissue, "up", negate = TRUE)
  FL_G_Tissue <- str_subset(FL_G_Tissue, "down", negate = TRUE)
  
  FL_G_Tissue <- str_sort(FL_G_Tissue, numeric = TRUE)
}
#clip off the .1 off the gene ids
de_genes_zf_ze_vs_fmg <- sapply(1:length(FL_G_Tissue),function(i){
  
  GeneL <- read.csv(FL_G_Tissue[i])$X
  GeneL <- as.character(GeneL)
  GeneL <- str_replace(GeneL,"[.]1","")
  
  print(GeneL)
  return(GeneL)
  
})
#compare them to TF baselist
de_genes_zf_ze_vs_fmg_TFs <- sapply(1:length(FL_G_Tissue),function(i){
  print(str_c("B",i))
  x <- intersect(unlist(de_genes_zf_ze_vs_fmg[i]),TF_Baselist)
  print(x)
  return(x)
})

view(de_genes_zf_fmg_vs_time_TFs)

Reduce(cbind,de_genes_zf_fmg_vs_time_TFs)

de_genes_zf_fmg_vs_time_TFs
de_genes_zf_ze_vs_time_TFs
de_genes_zf_ze_vs_fmg_TFs

write.csv(de_genes_zf_fmg_vs_time_TFs, file=here("analysis/TFs/de_genes_zf_fmg_vs_time_TFs.csv"))
write.csv(de_genes_zf_ze_vs_time_TFs, file=here("analysis/TFs/de_genes_zf_ze_vs_time_TFs.csv"))
write.csv(de_genes_zf_ze_vs_fmg_TFs, file=here("analysis/TFs/de_genes_zf_ze_vs_fmg_TFs.csv"))

for(i in 1:9){
  nam <- str_c("analysis/TFs/de_genes_zf_fmg_vs_time_B",i+1,"_vs_B",i,"_TFs.csv")
  de_genes_zf_fmg_vs_time_TFs[i]
  write.csv(de_genes_zf_fmg_vs_time_TFs[i], file=here(nam))
}

for(i in 1:9){
  nam <- str_c("analysis/TFs/de_genes_zf_ze_vs_time_B",i+1,"_vs_B",i,"_TFs.csv")
  de_genes_zf_ze_vs_time_TFs[i]
  write.csv(de_genes_zf_ze_vs_time_TFs[i], file=here(nam))
}

for(i in 1:10){
  nam <- str_c("analysis/TFs/de_genes_zf_ze_vs_fmg_B",i,"_TFs.csv")
  write.csv(de_genes_zf_ze_vs_time_TFs[i], file=here(nam))
}


intersect(topClusters$Cluster1,tg1_res) #12 genes


names(topClusters)
for(i in 1:length(topClusters)){
#  print(topClusters[i])
  print(names(topClusters[i]))
  cnam <- names(topClusters[i])
#  print(topClusters[i])
  print(intersect(topClusters[i],tg1_res))
#  print(intersect(topClusters$cnam,tg1_res))
}

topClusters_PC1_TFs

for(i in 1:length(topClusters)){
  #  print(topClusters[i])
#  print(names(topClusters[i]))
  cnam <- names(topClusters[i])
  #  print(topClusters[i])
  print(cnam)
  print(intersect(unlist(topClusters[i]),tg1_res))
  #  print(intersect(topClusters$cnam,tg1_res))
}


#loop above wont want to loop through and give me results as below, either in $ or [i] variants
#clusters below are the only ones that gave me genes
{tg1_res #81 total genes
intersect(topClusters$Cluster1,tg1_res) #12 genes
intersect(topClusters$Cluster5,tg1_res) #5 genes
intersect(topClusters$Cluster6,tg1_res) #2 genes
intersect(topClusters$Cluster7,tg1_res) #2 genes
intersect(topClusters$Cluster8,tg1_res) #1 gene
intersect(topClusters$Cluster9,tg1_res) #3 genes
intersect(topClusters$Cluster10,tg1_res) #1 gene
intersect(topClusters$Cluster11,tg1_res) #1 gene
intersect(topClusters$Cluster12,tg1_res) #12 genes
intersect(topClusters$Cluster14,tg1_res) #3 genes
intersect(topClusters$Cluster18,tg1_res) #1 gene
intersect(topClusters$Cluster19,tg1_res) #4 genes
intersect(topClusters$Cluster22,tg1_res) #1 gene
intersect(topClusters$Cluster23,tg1_res) #2 genes
intersect(topClusters$Cluster24,tg1_res) #1 gene
intersect(topClusters$Cluster25,tg1_res) #10 genes
intersect(topClusters$Cluster26,tg1_res) #2 genes
intersect(topClusters$Cluster28,tg1_res) #8 genes
intersect(topClusters$Cluster32,tg1_res) #1 gene
}
#disregard, I have fixed it


#export the TF genes from PC1, PC2 and PC3.



#FDN of TFs
source(here("Rtoolbox/src/getFDN.R"))

#READ FILE
geneList <- read_csv(here("analysis/DE/FMG_Time_B2_vs_B(-1)genes.csv"))
geneList <- geneList$X1
geneList <- str_replace_all(geneList, "[.]1","")

#READ SEIDR EDGELIST
edgeList <- read.table(here("data/seidr/network/cluster/edgeList.txt"),header=F,sep="\t")
edgeList <- edgeList[,1:2]

# data setup
# edgelist (2 columns, source and linked gene - 3rd column with value is removed)
# genes 


#subscript out of bounds error - why is that? Trying to access something that is out of the range of the matrix
#the error was that it was in a matrix


de_genes_zf_fmg_vs_time_TFs_FDN <- sapply(1:length(de_genes_zf_fmg_vs_time_TFs),function(i){
  print(str_c("B",i))
#  x <- intersect(unlist(de_genes_zf_ze_vs_time[i]),TF_Baselist)
  if(length(de_genes_zf_fmg_vs_time_TFs[[i]]) == 0){
    print("no genes")
    return(NULL)
  }else{
    x <- getFDN(edgeList,unlist(de_genes_zf_fmg_vs_time_TFs[i]))
    print(x)
    return(x)
  }
})

#FDNs
de_genes_zf_fmg_vs_time_TFs_FDN[[1]][,1]
de_genes_zf_fmg_vs_time_TFs_Transition <- sapply(1:length(de_genes_zf_fmg_vs_time_TFs),function(i){
  if(i == length(de_genes_zf_fmg_vs_time_TFs)){
    print("no further comparison available")
  }else{
    x <- intersect(unlist(de_genes_zf_fmg_vs_time_TFs_FDN[[i]][,1]),unlist(de_genes_zf_fmg_vs_time_TFs[i+1]))
    print(str_c("B",i+1," vs B",i, " intersecting FDN with B",i+2, " vs B", i+1," TFs"))
    print(x)
  }
})
de_genes_zf_fmg_vs_time_TFs_Transition_with_FDN <- sapply(1:length(de_genes_zf_fmg_vs_time_TFs),function(i){
  if(i == length(de_genes_zf_fmg_vs_time_TFs)){
    print("no further comparison available")
  }else{
  x <- intersect(unlist(de_genes_zf_fmg_vs_time_TFs_FDN[[i]][,1]),unlist(de_genes_zf_fmg_vs_time_TFs_FDN[[i+1]][,1]))
  print(str_c("B",i+1," vs B",i, " intersecting FDN with B",i+2, " vs B", i+1," TFs"))
  print(x)
  }
})
#WE HAVE TFS THAT TRANSITION FROM ONE TIME POINT TO ANOTHER

de_genes_zf_ze_vs_time_TFs_FDN <- sapply(1:length(de_genes_zf_ze_vs_time_TFs),function(i){
  print(str_c("B",i))
  #  x <- intersect(unlist(de_genes_zf_ze_vs_time[i]),TF_Baselist)
  if(length(de_genes_zf_ze_vs_time_TFs[[i]]) == 0){
    print("no genes")
    return(NULL)
  }else{
    x <- getFDN(edgeList,unlist(de_genes_zf_ze_vs_time_TFs[i]))
    print(x)
    return(x)
  }
})
#FDNs
de_genes_zf_ze_vs_time_TFs_Transition <- sapply(1:length(de_genes_zf_ze_vs_time_TFs),function(i){
  
  x <- intersect(unlist(de_genes_zf_ze_vs_time_TFs_FDN[[i]][,1]),unlist(de_genes_zf_ze_vs_time_TFs[i+1]))
  print(str_c("B",i," to B",i+1," TFs"))
  print(x)
  
})
#WE HAVE TFS THAT TRANSITION FROM ONE TIME POINT TO ANOTHER
de_genes_zf_ze_vs_time_TFs_Transition_with_FDN <- sapply(1:length(de_genes_zf_ze_vs_time_TFs),function(i){
  if(i == length(de_genes_zf_ze_vs_time_TFs)){
    print("no further comparison available")
  }else{
    x <- intersect(unlist(de_genes_zf_ze_vs_time_TFs_FDN[[i]][,1]),unlist(de_genes_zf_ze_vs_time_TFs_FDN[[i+1]][,1]))
    print(str_c("B",i+1," vs B",i, " intersecting FDN with B",i+2, " vs B", i+1," TFs"))
    print(x)
  }
})




#ZE VS FMG TISSUE TYPES
de_genes_zf_ze_vs_fmg_TFs_FDN <- sapply(1:length(de_genes_zf_ze_vs_fmg_TFs),function(i){
  print(str_c("B",i))
  #  x <- intersect(unlist(de_genes_zf_ze_vs_fmg[i]),TF_Baselist)
  if(length(de_genes_zf_ze_vs_fmg_TFs[[i]]) == 0){
    print("no genes")
    return(NULL)
  }else{
    x <- getFDN(edgeList,unlist(de_genes_zf_ze_vs_fmg_TFs[i]))
    print(x)
    return(x)
  }
})
#FDNs
de_genes_zf_ze_vs_fmg_TFs_Transition <- sapply(1:length(de_genes_zf_ze_vs_fmg_TFs),function(i){
  
  x <- intersect(unlist(de_genes_zf_ze_vs_fmg_TFs_FDN[[i]][,1]),unlist(de_genes_zf_ze_vs_fmg_TFs[i+1]))
  print(str_c("B",i," to B",i+1," TFs"))
  print(x)
  
})
#WE HAVE TFS THAT TRANSITION FROM ONE fmg POINT TO ANOTHER
de_genes_zf_ze_vs_fmg_TFs_Transition_with_FDN <- sapply(1:length(de_genes_zf_ze_vs_fmg_TFs),function(i){
  if(i == length(de_genes_zf_ze_vs_fmg_TFs)){
    print("no further comparison available")
  }else{
    x <- intersect(unlist(de_genes_zf_ze_vs_fmg_TFs_FDN[[i]][,1]),unlist(de_genes_zf_ze_vs_fmg_TFs_FDN[[i+1]][,1]))
    print(str_c("B",i+1," vs B",i, " intersecting FDN with B",i+2, " vs B", i+1," TFs"))
    print(x)
  }
})










de_genes_zf_fmg_vs_time_TFs_Transition
de_genes_zf_ze_vs_time_TFs_Transition

#ZF TIME FOR FMG
#ZF TIME FOR ZE
#ZF TISSUE EFFECT FOLDER

#SAME FOR SZ TIME FOR SE
#SAME FOR SZ TIME FOR ZE
#SAME FOR SZ TISSUE EFFECT FOLDER































#Compare if TFs are present in this DE list
#again, loop through, using intersect with TF_Baselist

#Extract only TFs that match in both

#Export file containing DE information of those TFs



### Do this for DE files as well as genes that contribute the most to the different Principle Components