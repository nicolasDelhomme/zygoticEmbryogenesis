library(igraph)



#READ SEIDR EDGELIST
edgeList <- read.table(here("data/seidr/network/cluster/edgeList.txt"),header=F,sep="\t")
edgeList <- edgeList[,1:2]


de_genes_zf_fmg_vs_time_TFs_Transition_with_FDN[1]

subset(edgeList, edgeList[,1] == unlist(de_genes_zf_fmg_vs_time_TFs_Transition_with_FDN[1]))


unlist(de_genes_zf_fmg_vs_time_TFs_Transition_with_FDN[1])


set.vertex.attribute(g1, "Origin", value = Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[1]][,1])
set.edge.attribute()

#FMG vs TIME
Edgelist_FDN_genes_zf_fmg_vs_time_TFs <- lapply(1:length(de_genes_zf_fmg_vs_time_TFs_Transition_with_FDN),function(i){
  
  return(edgeList %>% filter(V1 %in% unlist(de_genes_zf_fmg_vs_time_TFs_Transition_with_FDN[i])))
})

#FMG vs Time TFs with FDN network assemblies
for(i in 1:length(Edgelist_FDN_genes_zf_fmg_vs_time_TFs)){
  if(nrow(Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[i]]) == 0){
    print("no genes")
  }else{
    nam <- str_c("zf_fmg_vs_time_B",i+1,"B",i,"_vs_B",i+2,"B",i+1,"_TFs_with_FDNs_network")
    
    g <- graph.edgelist(as.matrix(Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[i]]), directed = FALSE)
    dir <- str_c("data/seidr/network/",nam)
    print(dir)
    write.graph(g, here(str_c("data/seidr/network/",nam)), format = c("ncol"))
  }
  if(i == length(Edgelist_FDN_genes_zf_fmg_vs_time_TFs)){
    nam <- str_c("zf_fmg_vs_time_TF_FDNs_unionifed_network")
    print(nam)
    inter <- bind_rows(Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[1]],
                       Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[2]],
                       Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[3]],
                       Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[4]],
                       Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[5]],
                       Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[6]],
                       Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[7]],
                       Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[8]],
                       Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[9]])
    inter <- distinct(inter)
    inter <- as.matrix(inter)
    inter <- graph.edgelist(inter, directed = FALSE)
    write.graph(inter, here(str_c("data/seidr/network/",nam)), format = c("ncol"))
  }
}


inter <- bind_rows(Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[1]],
                   Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[2]],
                   Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[3]],
                   Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[4]],
                   Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[5]],
                   Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[6]],
                   Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[7]],
                   Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[8]],
                   Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[9]])
l1 <- union(inter[,1],inter[,2])
l1 <- cbind(l1, rep("FMG_TIME",length(l1)))
rownames(l1) <- l1[,1]
l1 <- l1[,2]

write.csv(l1, file=here("data/seidr/network/zf_fmg_vs_time_TF_FDNs_unionifed_list.tsv"))

#ZE vs TIME
Edgelist_FDN_genes_zf_ze_vs_time_TFs <- lapply(1:length(de_genes_zf_ze_vs_time_TFs_Transition_with_FDN),function(i){
  
  return(edgeList %>% filter(V1 %in% unlist(de_genes_zf_ze_vs_time_TFs_Transition_with_FDN[i])))
})

Edgelist_FDN_genes_zf_ze_vs_time_TFs[[1]]


graphtest <- graph.edgelist(as.matrix(Edgelist_FDN_genes_zf_ze_vs_time_TFs[[5]]), directed = FALSE)
V(graphtest)

#ZE vs Time TFs with FDN network assemblies
for(i in 1:length(Edgelist_FDN_genes_zf_ze_vs_time_TFs)){
  if(nrow(Edgelist_FDN_genes_zf_ze_vs_time_TFs[[i]]) == 0){
    print("no genes")
  }else{
    nam <- str_c("zf_ze_vs_time_B",i+1,"B",i,"_vs_B",i+2,"B",i+1,"_TFs_with_FDNs_network")
    
    g <- graph.edgelist(as.matrix(Edgelist_FDN_genes_zf_ze_vs_time_TFs[[i]]), directed = FALSE)
    dir <- str_c("data/seidr/network/",nam)
    print(dir)
    write.graph(g, here(str_c("data/seidr/network/",nam)), format = c("ncol"))
  }
  if(i == length(Edgelist_FDN_genes_zf_ze_vs_time_TFs)){
    nam <- str_c("zf_ze_vs_time_TF_FDNs_unionifed_network")
    print(nam)
    inter <- bind_rows(Edgelist_FDN_genes_zf_ze_vs_time_TFs[[1]],
                       Edgelist_FDN_genes_zf_ze_vs_time_TFs[[2]],
                       Edgelist_FDN_genes_zf_ze_vs_time_TFs[[3]],
                       Edgelist_FDN_genes_zf_ze_vs_time_TFs[[4]],
                       Edgelist_FDN_genes_zf_ze_vs_time_TFs[[5]],
                       Edgelist_FDN_genes_zf_ze_vs_time_TFs[[6]],
                       Edgelist_FDN_genes_zf_ze_vs_time_TFs[[7]],
                       Edgelist_FDN_genes_zf_ze_vs_time_TFs[[8]],
                       Edgelist_FDN_genes_zf_ze_vs_time_TFs[[9]])
    inter <- distinct(inter)
    inter <- as.matrix(inter)
    inter <- graph.edgelist(inter, directed = FALSE)
    write.graph(inter, here(str_c("data/seidr/network/",nam)), format = c("ncol"))
  }
}


inter <- bind_rows(Edgelist_FDN_genes_zf_ze_vs_time_TFs[[1]],
                   Edgelist_FDN_genes_zf_ze_vs_time_TFs[[2]],
                   Edgelist_FDN_genes_zf_ze_vs_time_TFs[[3]],
                   Edgelist_FDN_genes_zf_ze_vs_time_TFs[[4]],
                   Edgelist_FDN_genes_zf_ze_vs_time_TFs[[5]],
                   Edgelist_FDN_genes_zf_ze_vs_time_TFs[[6]],
                   Edgelist_FDN_genes_zf_ze_vs_time_TFs[[7]],
                   Edgelist_FDN_genes_zf_ze_vs_time_TFs[[8]],
                   Edgelist_FDN_genes_zf_ze_vs_time_TFs[[9]])
l2 <- union(inter[,1],inter[,2])
l2 <- cbind(l2, rep("ZE_TIME",length(l2)))
rownames(l2) <- l2[,1]
l2 <- l2[,2]
write.csv(l2, file=here("data/seidr/network/zf_ze_vs_time_TF_FDNs_unionifed_list.tsv"))

#ZE vs FMG
Edgelist_FDN_genes_zf_ze_vs_fmg_TFs <- lapply(1:length(de_genes_zf_ze_vs_fmg_TFs_Transition_with_FDN),function(i){
  
  return(edgeList %>% filter(V1 %in% unlist(de_genes_zf_ze_vs_fmg_TFs_Transition_with_FDN[i])))
})

#ZE vs Time TFs with FDN network assemblies
for(i in 1:length(Edgelist_FDN_genes_zf_ze_vs_fmg_TFs)){
  if(nrow(Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[i]]) == 0){
    print("no genes")
  }else{
    nam <- str_c("zf_ze_vs_fmg_B",i+1,"B",i,"_vs_B",i+2,"B",i+1,"_TFs_with_FDNs_network")
    
    g <- graph.edgelist(as.matrix(Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[i]]), directed = FALSE)
    dir <- str_c("data/seidr/network/",nam)
    print(dir)
    write.graph(g, here(str_c("data/seidr/network/",nam)), format = c("ncol"))
  }
  if(i == length(Edgelist_FDN_genes_zf_ze_vs_fmg_TFs)){
    nam <- str_c("zf_ze_vs_fmg_TF_FDNs_unionifed_network")
    print(nam)
    inter <- bind_rows(Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[1]],
                       Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[2]],
                       Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[3]],
                       Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[4]],
                       Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[5]],
                       Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[6]],
                       Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[7]],
                       Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[8]],
                       Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[9]])
    inter <- distinct(inter)
    inter <- as.matrix(inter)
    inter <- graph.edgelist(inter, directed = FALSE)
    write.graph(inter, here(str_c("data/seidr/network/",nam)), format = c("ncol"))
  }
}


inter <- bind_rows(Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[1]],
                   Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[2]],
                   Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[3]],
                   Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[4]],
                   Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[5]],
                   Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[6]],
                   Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[7]],
                   Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[8]],
                   Edgelist_FDN_genes_zf_ze_vs_fmg_TFs[[9]])
l3 <- union(inter[,1],inter[,2])
l3 <- cbind(l3, rep("TISSUE_ZE_FMG",length(l3)))
rownames(l3) <- l3[,1]
l3 <- l3[,2]
write.csv(l3, file=here("data/seidr/network/zf_ze_vs_fmg_TF_FDNs_unionifed_list.tsv"))





l1.1 <- data.frame(names(l1),l1)
colnames(l1.1) <- c("GENE","FMG_TIME")
rownames(l1.1) <- 1:nrow(l1.1)
l2.1 <- data.frame(names(l2),l2)
colnames(l2.1) <- c("GENE","ZE_TIME")
rownames(l2.1) <- 1:nrow(l2.1)
l3.1 <- data.frame(names(l3),l3)
colnames(l3.1) <- c("GENE","ZE_FMG")
rownames(l3.1) <- 1:nrow(l3.1)

l4 <- full_join(l1.1,l2.1)
l4 <- full_join(l4,l3.1)
agg <- str_c(str_replace_na(l4[,2])," ",str_replace_na(l4[,3])," ",str_replace_na(l4[,4]))
agg <- str_replace(agg,"NA", "MISSING")


l5 <- cbind(l4[,1],agg)
l5
write.csv(l5, file=here("data/seidr/network/zf_ze_vs_fmg_TF_FDNs_unionifed_list_ALLMERGED_v2.tsv"))
#MA_10066964g0010



l1.1




rownames(l1.1) <- 1:nrow(l1.1)
l2.1 <- cbind(names(l2),l2)
rownames(l2.1) <- 1:nrow(l2.1)
l3.1 <- cbind(names(l3),l3)
rownames(l3.1) <- 1:nrow(l3.1)

l1.1 <- as.data.frame(l1.1)
l2.1 <- as.data.frame(l2.1)
l3.1 <- as.data.frame(l3.1)

l4 <- full_join(l1.1,l2.1)
l4 <- full_join(l4,l3.1)
l4



Edgelist_FDN_genes_zf_fmg_vs_time_TFs[[1]]

ed_mat <- edgeList %>% filter(V1 %in% unlist(de_genes_zf_fmg_vs_time_TFs_Transition_with_FDN[1]))
ed_mat <- as.matrix(ed_mat)
class(ed_mat)


graph.edgelist(ed_mat, directed = FALSE)


