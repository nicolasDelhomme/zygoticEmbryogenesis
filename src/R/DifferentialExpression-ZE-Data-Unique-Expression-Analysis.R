#Differential Gene Expression Identification of Unique Genes under different conditions in Zygotic Embryogenesis Dataset

FL_FMG_Time <- list.files(here("analysis/DE/ZF_FMG_Time/"), 
                          recursive = TRUE, 
                          pattern = "FMG_Time",
                          full.names = TRUE)

FL_FMG_Time <- str_subset(FL_FMG_Time, "genes.csv")
FL_FMG_Time <- str_subset(FL_FMG_Time, "up", negate = TRUE)
FL_FMG_Time <- str_subset(FL_FMG_Time, "down", negate = TRUE)

FMG_GL <- sapply(1:length(FL_FMG_Time), function(i){

  a <- read.csv(FL_FMG_Time[i])$X
  a <- as.character(a)
  a <- str_replace(a,"[.]1","")
  return(a)
})

length(unlist(FMG_GL))

rownames(unlist(FMG_GL))
unique(unlist(FMG_GL))

FMG_GL_Names <- sapply(1:length(FL_FMG_Time), function(i){
  a <- str_split(FL_FMG_Time, "FMG_")[[i]][3]
  a <- str_replace(a, ".csv","")
  return(a)
})
#FMG_GL_Names <- str_sort(FMG_GL_Names, numeric = TRUE)
names(FMG_GL) <- FMG_GL_Names

FMG_GL_Uniq <- sapply(1:length(FMG_GL),function(i){
  if(i == 1){
    print((i+1):length(FMG_GL))
    a <- ((i+1):length(FMG_GL))
    
    print("BREAK")
    reference <- unlist(FMG_GL[a])
    
    #get uniq
    uniq <- FMG_GL[i]
    uniq <- setdiff(uniq,reference)
    names(uniq) <- names(FMG_GL)[i]
    return(uniq)
    
    
  }else if(i == 9){
    print(1:(i-1))
    a <- (1:(i-1))
    
    print("BREAK")
    reference <- unlist(FMG_GL[a])
    
    #get uniq
    uniq <- FMG_GL[i]
    uniq <- setdiff(uniq,reference)
    names(uniq) <- names(FMG_GL)[i]
    return(uniq)
    
  }else{
    print(1:(i-1))
    a <- (1:(i-1))
    
    print((i+1):length(FMG_GL))
    b <- ((i+1):length(FMG_GL))
    
    print("BREAK")  
    reference <- unlist(c(FMG_GL[a],FMG_GL[b]))
    
    #get uniq
    uniq <- FMG_GL[i]
    uniq <- setdiff(uniq,reference)
    names(uniq) <- names(FMG_GL)[i]
    return(uniq)
      }
})



#looking at Time in FMG Tissue
#B2vsB1
FMG_B2_vs_B1_GL <- read.csv(FL_FMG_Time[2])$X
FMG_B2_vs_B1_GL <- FMG_B2_vs_B1_GL$X
FMG_B2_vs_B1_GL <- as.character(FMG_B2_vs_B1_GL)
#B3vsB2
FMG_B3_vs_B2_GL <- read.csv(FL_FMG_Time[3])
FMG_B3_vs_B2_GL <- FMG_B3_vs_B2_GL$X
FMG_B3_vs_B2_GL <- as.character(FMG_B3_vs_B2_GL)
#B4vsB3
FMG_B4_vs_B3_GL <- read.csv(FL_FMG_Time[4])
FMG_B4_vs_B3_GL <- FMG_B4_vs_B3_GL$X
FMG_B4_vs_B3_GL <- as.character(FMG_B4_vs_B3_GL)
#B5vsB4
FMG_B5_vs_B4_GL <- read.csv(FL_FMG_Time[5])
FMG_B5_vs_B4_GL <- FMG_B5_vs_B4_GL$X
FMG_B5_vs_B4_GL <- as.character(FMG_B5_vs_B4_GL)
#B6vsB5
FMG_B6_vs_B5_GL <- read.csv(FL_FMG_Time[6])
FMG_B6_vs_B5_GL <- FMG_B6_vs_B5_GL$X
FMG_B6_vs_B5_GL <- as.character(FMG_B6_vs_B5_GL)
#B7vsB6
FMG_B7_vs_B6_GL <- read.csv(FL_FMG_Time[7])
FMG_B7_vs_B6_GL <- FMG_B7_vs_B6_GL$X
FMG_B7_vs_B6_GL <- as.character(FMG_B7_vs_B6_GL)
#B8vsB7
FMG_B8_vs_B7_GL <- read.csv(FL_FMG_Time[8])
FMG_B8_vs_B7_GL <- FMG_B8_vs_B7_GL$X
FMG_B8_vs_B7_GL <- as.character(FMG_B8_vs_B7_GL)
#B9vsB8
FMG_B9_vs_B8_GL <- read.csv(FL_FMG_Time[9])
FMG_B9_vs_B8_GL <- FMG_B9_vs_B8_GL$X
FMG_B9_vs_B8_GL <- as.character(FMG_B9_vs_B8_GL)
#B10vsB9
FMG_B10_vs_B9_GL <- read.csv(FL_FMG_Time[1])
FMG_B10_vs_B9_GL <- FMG_B10_vs_B9_GL$X
FMG_B10_vs_B9_GL <- as.character(FMG_B10_vs_B9_GL)

###compare every GL and keep only the unique ones
c(FMG_B2_vs_B1_GL,FMG_B3_vs_B2_GL)

FMG_GL_Repository$Total <- list(FMG_B2_vs_B1_GL,FMG_B3_vs_B2_GL,FMG_B4_vs_B3_GL,
                          FMG_B5_vs_B4_GL,FMG_B6_vs_B5_GL,FMG_B7_vs_B6_GL,
                          FMG_B8_vs_B7_GL,FMG_B9_vs_B8_GL,FMG_B10_vs_B9_GL)


#combine all GLs except for the one in question. 

REFERENCE <- unlist(FMG_GL_Repository[2:9])
FMG_B2_vs_B1_GL_Uniq <- setdiff(FMG_B2_vs_B1_GL,REFERENCE)
FMG_GL_Repository$B2 <- setdiff(FMG_B2_vs_B1_GL,REFERENCE)

REFERENCE <- unlist(c(FMG_GL_Repository[1],FMG_GL_Repository[3:9]))
FMG_GL_Repository$B3 <- setdiff(FMG_B3_vs_B2_GL,REFERENCE)

REFERENCE <- unlist(c(FMG_GL_Repository[1:2],FMG_GL_Repository[4:9]))
FMG_GL_Repository$B4 <- setdiff(FMG_B4_vs_B3_GL,REFERENCE)

REFERENCE <- unlist(c(FMG_GL_Repository[1:3],FMG_GL_Repository[5:9]))
FMG_GL_Repository$B5 <- setdiff(FMG_B5_vs_B4_GL,REFERENCE)

REFERENCE <- unlist(c(FMG_GL_Repository[1:4],FMG_GL_Repository[6:9]))
FMG_GL_Repository$B6 <- setdiff(FMG_B6_vs_B5_GL,REFERENCE)

REFERENCE <- unlist(c(FMG_GL_Repository[1:5],FMG_GL_Repository[7:9]))
FMG_GL_Repository$B7 <- setdiff(FMG_B7_vs_B6_GL,REFERENCE)

REFERENCE <- unlist(c(FMG_GL_Repository[1:6],FMG_GL_Repository[8:9]))
FMG_GL_Repository$B8 <- setdiff(FMG_B8_vs_B7_GL,REFERENCE)

REFERENCE <- unlist(c(FMG_GL_Repository[1:7],FMG_GL_Repository[9]))
FMG_GL_Repository$B9 <- setdiff(FMG_B9_vs_B8_GL,REFERENCE)

REFERENCE <- unlist(c(FMG_GL_Repository[1:8]))
FMG_GL_Repository$B10 <- setdiff(FMG_B10_vs_B9_GL,REFERENCE)


FMG_GL_Repository
length(FMG_GL_Repository)

function(i){
  REFERENCE <- unlist(c(FMG_GL_Repository[i:length(FMG_GL_Repository)]))
  FMG_GL_Repository$B10 <- setdiff(FMG_B10_vs_B9_GL,REFERENCE)
  
}





# gopher
source(here("UPSCb-common/src/R/gopher.R"))
background <- rownames(vst)
background <- str_replace(background,"[.]1","")
FMG_GL_Repository$B2 <- str_replace(FMG_GL_Repository$B2,"[.]1","")
FMG_GL_Repository$B3 <- str_replace(FMG_GL_Repository$B3,"[.]1","")
FMG_GL_Repository$B4 <- str_replace(FMG_GL_Repository$B4,"[.]1","")
FMG_GL_Repository$B5 <- str_replace(FMG_GL_Repository$B5,"[.]1","")
FMG_GL_Repository$B6 <- str_replace(FMG_GL_Repository$B6,"[.]1","")
FMG_GL_Repository$B7 <- str_replace(FMG_GL_Repository$B7,"[.]1","")
FMG_GL_Repository$B8 <- str_replace(FMG_GL_Repository$B8,"[.]1","")
FMG_GL_Repository$B9 <- str_replace(FMG_GL_Repository$B9,"[.]1","")
FMG_GL_Repository$B10 <- str_replace(FMG_GL_Repository$B10,"[.]1","")



FMG_GL_Repository$B2
GL_FMG_B2B1_enr <- gopher(genes = FMG_GL_Repository$B2, background = background,task = c("go","mapman"),url = "pabies")
GL_FMG_B3B2_enr <- gopher(genes = FMG_GL_Repository$B3, background = background,task = c("go","mapman"),url = "pabies")
GL_FMG_B4B3_enr <- gopher(genes = FMG_GL_Repository$B4, background = background,task = c("go","mapman"),url = "pabies")
GL_FMG_B5B4_enr <- gopher(genes = FMG_GL_Repository$B5, background = background,task = c("go","mapman"),url = "pabies")
GL_FMG_B6B5_enr <- gopher(genes = FMG_GL_Repository$B6, background = background,task = c("go","mapman"),url = "pabies")
GL_FMG_B7B6_enr <- gopher(genes = FMG_GL_Repository$B7, background = background,task = c("go","mapman"),url = "pabies")
GL_FMG_B8B7_enr <- gopher(genes = FMG_GL_Repository$B8, background = background,task = c("go","mapman"),url = "pabies")
GL_FMG_B9B8_enr <- gopher(genes = FMG_GL_Repository$B9, background = background,task = c("go","mapman"),url = "pabies")
GL_FMG_B10B9_enr <- gopher(genes = FMG_GL_Repository$B10, background = background,task = c("go","mapman"),url = "pabies")


GL_FMG_B2B1_enr
GL_FMG_B3B2_enr
GL_FMG_B4B3_enr
GL_FMG_B5B4_enr
GL_FMG_B6B5_enr
GL_FMG_B7B6_enr
GL_FMG_B8B7_enr
GL_FMG_B9B8_enr
GL_FMG_B10B9_enr



source(here("Rtoolbox/src/plotEnrichedTreemap.R"))

#FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_FMG_B2B1_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_FMG_B2B1_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_FMG_B3B2_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_FMG_B3B2_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_FMG_B5B4_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_FMG_B5B4_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_FMG_B6B5_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_FMG_B6B5_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_FMG_B7B6_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_FMG_B7B6_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_FMG_B8B7_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_FMG_B8B7_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_FMG_B9B8_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_FMG_B9B8_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_FMG_B10B9_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_FMG_B10B9_enr, enrichment = "mapman", clusterColor = "#1B75BC")



#ZE vs Time
FL_ZE_Time <- list.files(here("analysis/DE"), 
                          recursive = TRUE, 
                          pattern = "ZE_Time",
                          full.names = TRUE)

FL_ZE_Time <- str_subset(FL_ZE_Time, "genes.csv")
FL_ZE_Time <- str_subset(FL_ZE_Time, "up", negate = TRUE)
FL_ZE_Time <- str_subset(FL_ZE_Time, "down", negate = TRUE)

#looking at Time in ZE Tissue
#B2vsB1
ZE_B2_vs_B1_GL <- read.csv(FL_ZE_Time[2])$X
ZE_B2_vs_B1_GL <- ZE_B2_vs_B1_GL$X
ZE_B2_vs_B1_GL <- as.character(ZE_B2_vs_B1_GL)
#B3vsB2
ZE_B3_vs_B2_GL <- read.csv(FL_ZE_Time[3])
ZE_B3_vs_B2_GL <- ZE_B3_vs_B2_GL$X
ZE_B3_vs_B2_GL <- as.character(ZE_B3_vs_B2_GL)
#B4vsB3
ZE_B4_vs_B3_GL <- read.csv(FL_ZE_Time[4])
ZE_B4_vs_B3_GL <- ZE_B4_vs_B3_GL$X
ZE_B4_vs_B3_GL <- as.character(ZE_B4_vs_B3_GL)
#B5vsB4
ZE_B5_vs_B4_GL <- read.csv(FL_ZE_Time[5])
ZE_B5_vs_B4_GL <- ZE_B5_vs_B4_GL$X
ZE_B5_vs_B4_GL <- as.character(ZE_B5_vs_B4_GL)
#B6vsB5
ZE_B6_vs_B5_GL <- read.csv(FL_ZE_Time[6])
ZE_B6_vs_B5_GL <- ZE_B6_vs_B5_GL$X
ZE_B6_vs_B5_GL <- as.character(ZE_B6_vs_B5_GL)
#B7vsB6
ZE_B7_vs_B6_GL <- read.csv(FL_ZE_Time[7])
ZE_B7_vs_B6_GL <- ZE_B7_vs_B6_GL$X
ZE_B7_vs_B6_GL <- as.character(ZE_B7_vs_B6_GL)
#B8vsB7
ZE_B8_vs_B7_GL <- read.csv(FL_ZE_Time[8])
ZE_B8_vs_B7_GL <- ZE_B8_vs_B7_GL$X
ZE_B8_vs_B7_GL <- as.character(ZE_B8_vs_B7_GL)
#B9vsB8
ZE_B9_vs_B8_GL <- read.csv(FL_ZE_Time[9])
ZE_B9_vs_B8_GL <- ZE_B9_vs_B8_GL$X
ZE_B9_vs_B8_GL <- as.character(ZE_B9_vs_B8_GL)
#B10vsB9
ZE_B10_vs_B9_GL <- read.csv(FL_ZE_Time[1])
ZE_B10_vs_B9_GL <- ZE_B10_vs_B9_GL$X
ZE_B10_vs_B9_GL <- as.character(ZE_B10_vs_B9_GL)

###compare every GL and keep only the unique ones
c(ZE_B2_vs_B1_GL,ZE_B3_vs_B2_GL)
ZE_GL_Repository <- list(NULL)
ZE_GL_Repository <- list(ZE_B2_vs_B1_GL,ZE_B3_vs_B2_GL,ZE_B4_vs_B3_GL,
                                ZE_B5_vs_B4_GL,ZE_B6_vs_B5_GL,ZE_B7_vs_B6_GL,
                                ZE_B8_vs_B7_GL,ZE_B9_vs_B8_GL,ZE_B10_vs_B9_GL)


#combine all GLs except for the one in question. 
{
REFERENCE <- unlist(ZE_GL_Repository[2:9])
ZE_B2_vs_B1_GL_Uniq <- setdiff(ZE_B2_vs_B1_GL,REFERENCE)
ZE_GL_Repository$B2 <- setdiff(ZE_B2_vs_B1_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_GL_Repository[1],ZE_GL_Repository[3:9]))
ZE_GL_Repository$B3 <- setdiff(ZE_B3_vs_B2_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_GL_Repository[1:2],ZE_GL_Repository[4:9]))
ZE_GL_Repository$B4 <- setdiff(ZE_B4_vs_B3_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_GL_Repository[1:3],ZE_GL_Repository[5:9]))
ZE_GL_Repository$B5 <- setdiff(ZE_B5_vs_B4_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_GL_Repository[1:4],ZE_GL_Repository[6:9]))
ZE_GL_Repository$B6 <- setdiff(ZE_B6_vs_B5_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_GL_Repository[1:5],ZE_GL_Repository[7:9]))
ZE_GL_Repository$B7 <- setdiff(ZE_B7_vs_B6_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_GL_Repository[1:6],ZE_GL_Repository[8:9]))
ZE_GL_Repository$B8 <- setdiff(ZE_B8_vs_B7_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_GL_Repository[1:7],ZE_GL_Repository[9]))
ZE_GL_Repository$B9 <- setdiff(ZE_B9_vs_B8_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_GL_Repository[1:8]))
ZE_GL_Repository$B10 <- setdiff(ZE_B10_vs_B9_GL,REFERENCE)
}

ZE_GL_Repository
length(ZE_GL_Repository)


# gopher
source(here("UPSCb-common/src/R/gopher.R"))
background <- rownames(vst)
#preprocessing
{background <- str_replace(background,"[.]1","")
ZE_GL_Repository$B2 <- str_replace(ZE_GL_Repository$B2,"[.]1","")
ZE_GL_Repository$B3 <- str_replace(ZE_GL_Repository$B3,"[.]1","")
ZE_GL_Repository$B4 <- str_replace(ZE_GL_Repository$B4,"[.]1","")
ZE_GL_Repository$B5 <- str_replace(ZE_GL_Repository$B5,"[.]1","")
ZE_GL_Repository$B6 <- str_replace(ZE_GL_Repository$B6,"[.]1","")
ZE_GL_Repository$B7 <- str_replace(ZE_GL_Repository$B7,"[.]1","")
ZE_GL_Repository$B8 <- str_replace(ZE_GL_Repository$B8,"[.]1","")
ZE_GL_Repository$B9 <- str_replace(ZE_GL_Repository$B9,"[.]1","")
ZE_GL_Repository$B10 <- str_replace(ZE_GL_Repository$B10,"[.]1","")
}

#enrichment of genes
{GL_ZE_B2B1_enr <- gopher(genes = ZE_GL_Repository$B2, background = background,task = c("go","mapman"),url = "pabies") #no genes
GL_ZE_B3B2_enr <- gopher(genes = ZE_GL_Repository$B3, background = background,task = c("go","mapman"),url = "pabies") #no genes
GL_ZE_B4B3_enr <- gopher(genes = ZE_GL_Repository$B4, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_B5B4_enr <- gopher(genes = ZE_GL_Repository$B5, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_B6B5_enr <- gopher(genes = ZE_GL_Repository$B6, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_B7B6_enr <- gopher(genes = ZE_GL_Repository$B7, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_B8B7_enr <- gopher(genes = ZE_GL_Repository$B8, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_B9B8_enr <- gopher(genes = ZE_GL_Repository$B9, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_B10B9_enr <- gopher(genes = ZE_GL_Repository$B10, background = background,task = c("go","mapman"),url = "pabies")
}

#check enrichment numbers
{GL_ZE_B2B1_enr
GL_ZE_B3B2_enr
GL_ZE_B4B3_enr
GL_ZE_B5B4_enr
GL_ZE_B6B5_enr
GL_ZE_B7B6_enr
GL_ZE_B8B7_enr
GL_ZE_B9B8_enr
GL_ZE_B10B9_enr
}

#PLotTreemaps
{#ZE_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_B2B1_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_B2B1_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_B3B2_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_B3B2_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_B5B4_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_B5B4_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_B6B5_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_B6B5_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_B7B6_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_B7B6_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_B8B7_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_B8B7_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_B9B8_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_B9B8_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_B10B9_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_B10B9_enr, enrichment = "mapman", clusterColor = "#1B75BC")
}





FL_ZE_vs_FMG_Tissue <- list.files(here("analysis/DE"), 
                          recursive = TRUE, 
                          pattern = "ZE_vs_FMG",
                          full.names = TRUE)

FL_ZE_vs_FMG_Tissue <- str_subset(FL_ZE_vs_FMG_Tissue, "genes.csv")
FL_ZE_vs_FMG_Tissue <- str_subset(FL_ZE_vs_FMG_Tissue, "up", negate = TRUE)
FL_ZE_vs_FMG_Tissue <- str_subset(FL_ZE_vs_FMG_Tissue, "down", negate = TRUE)

#looking at Tissue in ZE_vs_FMG Tissue
#B2vsB1
ZE_vs_FMG_B2_vs_B1_GL <- read.csv(FL_ZE_vs_FMG_Tissue[2])$X
ZE_vs_FMG_B2_vs_B1_GL <- ZE_vs_FMG_B2_vs_B1_GL$X
ZE_vs_FMG_B2_vs_B1_GL <- as.character(ZE_vs_FMG_B2_vs_B1_GL)
#B3vsB2
ZE_vs_FMG_B3_vs_B2_GL <- read.csv(FL_ZE_vs_FMG_Tissue[3])
ZE_vs_FMG_B3_vs_B2_GL <- ZE_vs_FMG_B3_vs_B2_GL$X
ZE_vs_FMG_B3_vs_B2_GL <- as.character(ZE_vs_FMG_B3_vs_B2_GL)
#B4vsB3
ZE_vs_FMG_B4_vs_B3_GL <- read.csv(FL_ZE_vs_FMG_Tissue[4])
ZE_vs_FMG_B4_vs_B3_GL <- ZE_vs_FMG_B4_vs_B3_GL$X
ZE_vs_FMG_B4_vs_B3_GL <- as.character(ZE_vs_FMG_B4_vs_B3_GL)
#B5vsB4
ZE_vs_FMG_B5_vs_B4_GL <- read.csv(FL_ZE_vs_FMG_Tissue[5])
ZE_vs_FMG_B5_vs_B4_GL <- ZE_vs_FMG_B5_vs_B4_GL$X
ZE_vs_FMG_B5_vs_B4_GL <- as.character(ZE_vs_FMG_B5_vs_B4_GL)
#B6vsB5
ZE_vs_FMG_B6_vs_B5_GL <- read.csv(FL_ZE_vs_FMG_Tissue[6])
ZE_vs_FMG_B6_vs_B5_GL <- ZE_vs_FMG_B6_vs_B5_GL$X
ZE_vs_FMG_B6_vs_B5_GL <- as.character(ZE_vs_FMG_B6_vs_B5_GL)
#B7vsB6
ZE_vs_FMG_B7_vs_B6_GL <- read.csv(FL_ZE_vs_FMG_Tissue[7])
ZE_vs_FMG_B7_vs_B6_GL <- ZE_vs_FMG_B7_vs_B6_GL$X
ZE_vs_FMG_B7_vs_B6_GL <- as.character(ZE_vs_FMG_B7_vs_B6_GL)
#B8vsB7
ZE_vs_FMG_B8_vs_B7_GL <- read.csv(FL_ZE_vs_FMG_Tissue[8])
ZE_vs_FMG_B8_vs_B7_GL <- ZE_vs_FMG_B8_vs_B7_GL$X
ZE_vs_FMG_B8_vs_B7_GL <- as.character(ZE_vs_FMG_B8_vs_B7_GL)
#B9vsB8
ZE_vs_FMG_B9_vs_B8_GL <- read.csv(FL_ZE_vs_FMG_Tissue[9])
ZE_vs_FMG_B9_vs_B8_GL <- ZE_vs_FMG_B9_vs_B8_GL$X
ZE_vs_FMG_B9_vs_B8_GL <- as.character(ZE_vs_FMG_B9_vs_B8_GL)
#B10vsB9
ZE_vs_FMG_B10_vs_B9_GL <- read.csv(FL_ZE_vs_FMG_Tissue[1])
ZE_vs_FMG_B10_vs_B9_GL <- ZE_vs_FMG_B10_vs_B9_GL$X
ZE_vs_FMG_B10_vs_B9_GL <- as.character(ZE_vs_FMG_B10_vs_B9_GL)

###compare every GL and keep only the unique ones
c(ZE_vs_FMG_B2_vs_B1_GL,ZE_vs_FMG_B3_vs_B2_GL)

ZE_vs_FMG_GL_Repository <- list(ZE_vs_FMG_B2_vs_B1_GL,ZE_vs_FMG_B3_vs_B2_GL,ZE_vs_FMG_B4_vs_B3_GL,
                                ZE_vs_FMG_B5_vs_B4_GL,ZE_vs_FMG_B6_vs_B5_GL,ZE_vs_FMG_B7_vs_B6_GL,
                                ZE_vs_FMG_B8_vs_B7_GL,ZE_vs_FMG_B9_vs_B8_GL,ZE_vs_FMG_B10_vs_B9_GL)


#combine all GLs except for the one in question. 
{
REFERENCE <- unlist(ZE_vs_FMG_GL_Repository[2:9])
ZE_vs_FMG_B2_vs_B1_GL_Uniq <- setdiff(ZE_vs_FMG_B2_vs_B1_GL,REFERENCE)
ZE_vs_FMG_GL_Repository$B2 <- setdiff(ZE_vs_FMG_B2_vs_B1_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_vs_FMG_GL_Repository[1],ZE_vs_FMG_GL_Repository[3:9]))
ZE_vs_FMG_GL_Repository$B3 <- setdiff(ZE_vs_FMG_B3_vs_B2_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_vs_FMG_GL_Repository[1:2],ZE_vs_FMG_GL_Repository[4:9]))
ZE_vs_FMG_GL_Repository$B4 <- setdiff(ZE_vs_FMG_B4_vs_B3_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_vs_FMG_GL_Repository[1:3],ZE_vs_FMG_GL_Repository[5:9]))
ZE_vs_FMG_GL_Repository$B5 <- setdiff(ZE_vs_FMG_B5_vs_B4_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_vs_FMG_GL_Repository[1:4],ZE_vs_FMG_GL_Repository[6:9]))
ZE_vs_FMG_GL_Repository$B6 <- setdiff(ZE_vs_FMG_B6_vs_B5_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_vs_FMG_GL_Repository[1:5],ZE_vs_FMG_GL_Repository[7:9]))
ZE_vs_FMG_GL_Repository$B7 <- setdiff(ZE_vs_FMG_B7_vs_B6_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_vs_FMG_GL_Repository[1:6],ZE_vs_FMG_GL_Repository[8:9]))
ZE_vs_FMG_GL_Repository$B8 <- setdiff(ZE_vs_FMG_B8_vs_B7_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_vs_FMG_GL_Repository[1:7],ZE_vs_FMG_GL_Repository[9]))
ZE_vs_FMG_GL_Repository$B9 <- setdiff(ZE_vs_FMG_B9_vs_B8_GL,REFERENCE)

REFERENCE <- unlist(c(ZE_vs_FMG_GL_Repository[1:8]))
ZE_vs_FMG_GL_Repository$B10 <- setdiff(ZE_vs_FMG_B10_vs_B9_GL,REFERENCE)
}

ZE_vs_FMG_GL_Repository
length(ZE_vs_FMG_GL_Repository)

function(i){
  REFERENCE <- unlist(c(ZE_vs_FMG_GL_Repository[i:length(ZE_vs_FMG_GL_Repository)]))
  ZE_vs_FMG_GL_Repository$B10 <- setdiff(ZE_vs_FMG_B10_vs_B9_GL,REFERENCE)
  
}





# gopher
source(here("UPSCb-common/src/R/gopher.R"))
background <- rownames(vst)
background <- str_replace(background,"[.]1","")
ZE_vs_FMG_GL_Repository$B2 <- str_replace(ZE_vs_FMG_GL_Repository$B2,"[.]1","")
ZE_vs_FMG_GL_Repository$B3 <- str_replace(ZE_vs_FMG_GL_Repository$B3,"[.]1","")
ZE_vs_FMG_GL_Repository$B4 <- str_replace(ZE_vs_FMG_GL_Repository$B4,"[.]1","")
ZE_vs_FMG_GL_Repository$B5 <- str_replace(ZE_vs_FMG_GL_Repository$B5,"[.]1","")
ZE_vs_FMG_GL_Repository$B6 <- str_replace(ZE_vs_FMG_GL_Repository$B6,"[.]1","")
ZE_vs_FMG_GL_Repository$B7 <- str_replace(ZE_vs_FMG_GL_Repository$B7,"[.]1","")
ZE_vs_FMG_GL_Repository$B8 <- str_replace(ZE_vs_FMG_GL_Repository$B8,"[.]1","")
ZE_vs_FMG_GL_Repository$B9 <- str_replace(ZE_vs_FMG_GL_Repository$B9,"[.]1","")
ZE_vs_FMG_GL_Repository$B10 <- str_replace(ZE_vs_FMG_GL_Repository$B10,"[.]1","")




ZE_vs_FMG_GL_Repository$B2
ZE_vs_FMG_GL_Repository$B3
ZE_vs_FMG_GL_Repository$B4
ZE_vs_FMG_GL_Repository$B5
ZE_vs_FMG_GL_Repository$B6
ZE_vs_FMG_GL_Repository$B7
ZE_vs_FMG_GL_Repository$B8
ZE_vs_FMG_GL_Repository$B9
ZE_vs_FMG_GL_Repository$B10


FMG_GL_Repository$B2
FMG_GL_Repository$B3
FMG_GL_Repository$B4
FMG_GL_Repository$B5
FMG_GL_Repository$B6
FMG_GL_Repository$B7
FMG_GL_Repository$B8
FMG_GL_Repository$B9
FMG_GL_Repository$B10

ZE_GL_Repository$B2
ZE_GL_Repository$B3
ZE_GL_Repository$B4
ZE_GL_Repository$B5
ZE_GL_Repository$B6
ZE_GL_Repository$B7
ZE_GL_Repository$B8
ZE_GL_Repository$B9
ZE_GL_Repository$B10






GL_ZE_vs_FMG_B2B1_enr <- gopher(genes = ZE_vs_FMG_GL_Repository$B2, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_vs_FMG_B3B2_enr <- gopher(genes = ZE_vs_FMG_GL_Repository$B3, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_vs_FMG_B4B3_enr <- gopher(genes = ZE_vs_FMG_GL_Repository$B4, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_vs_FMG_B5B4_enr <- gopher(genes = ZE_vs_FMG_GL_Repository$B5, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_vs_FMG_B6B5_enr <- gopher(genes = ZE_vs_FMG_GL_Repository$B6, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_vs_FMG_B7B6_enr <- gopher(genes = ZE_vs_FMG_GL_Repository$B7, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_vs_FMG_B8B7_enr <- gopher(genes = ZE_vs_FMG_GL_Repository$B8, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_vs_FMG_B9B8_enr <- gopher(genes = ZE_vs_FMG_GL_Repository$B9, background = background,task = c("go","mapman"),url = "pabies")
GL_ZE_vs_FMG_B10B9_enr <- gopher(genes = ZE_vs_FMG_GL_Repository$B10, background = background,task = c("go","mapman"),url = "pabies")



GL_ZE_vs_FMG_B2B1_enr
GL_ZE_vs_FMG_B3B2_enr
GL_ZE_vs_FMG_B4B3_enr
GL_ZE_vs_FMG_B5B4_enr
GL_ZE_vs_FMG_B6B5_enr
GL_ZE_vs_FMG_B7B6_enr
GL_ZE_vs_FMG_B8B7_enr
GL_ZE_vs_FMG_B9B8_enr
GL_ZE_vs_FMG_B10B9_enr

#ZE_vs_FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_vs_FMG_B2B1_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_vs_FMG_B2B1_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_vs_FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_vs_FMG_B3B2_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_vs_FMG_B3B2_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_vs_FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_vs_FMG_B5B4_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_vs_FMG_B5B4_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_vs_FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_vs_FMG_B6B5_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_vs_FMG_B6B5_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_vs_FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_vs_FMG_B7B6_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_vs_FMG_B7B6_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_vs_FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_vs_FMG_B8B7_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_vs_FMG_B8B7_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_vs_FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_vs_FMG_B9B8_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_vs_FMG_B9B8_enr, enrichment = "mapman", clusterColor = "#1B75BC")

#ZE_vs_FMG_B2 GO and Mapman Treemap
plotEnrichedTreemap(GL_ZE_vs_FMG_B10B9_enr, enrichment = "go", namespace = "none")
plotEnrichedTreemap(GL_ZE_vs_FMG_B10B9_enr, enrichment = "mapman", clusterColor = "#1B75BC")







#FMG TIME (UNIQUE)
GL_FMG_B2B1_enr
GL_FMG_B3B2_enr
GL_FMG_B4B3_enr
GL_FMG_B5B4_enr
GL_FMG_B6B5_enr
GL_FMG_B7B6_enr
GL_FMG_B8B7_enr
GL_FMG_B9B8_enr
GL_FMG_B10B9_enr

#PLOTS FMG TIME
{
  #FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_FMG_B2B1_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_FMG_B2B1_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_FMG_B3B2_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_FMG_B3B2_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_FMG_B5B4_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_FMG_B5B4_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_FMG_B6B5_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_FMG_B6B5_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_FMG_B7B6_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_FMG_B7B6_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_FMG_B8B7_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_FMG_B8B7_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_FMG_B9B8_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_FMG_B9B8_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_FMG_B10B9_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_FMG_B10B9_enr, enrichment = "mapman", clusterColor = "#1B75BC")
}


#ZE TIME (UNIQUE)
GL_ZE_B2B1_enr
GL_ZE_B3B2_enr
GL_ZE_B4B3_enr
GL_ZE_B5B4_enr
GL_ZE_B6B5_enr
GL_ZE_B7B6_enr
GL_ZE_B8B7_enr
GL_ZE_B9B8_enr
GL_ZE_B10B9_enr
 
#PLOTS ZE TIME
{
  #ZE_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_B2B1_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_B2B1_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_B3B2_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_B3B2_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_B5B4_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_B5B4_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_B6B5_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_B6B5_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_B7B6_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_B7B6_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_B8B7_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_B8B7_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_B9B8_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_B9B8_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_B10B9_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_B10B9_enr, enrichment = "mapman", clusterColor = "#1B75BC")
}

#ZE VS FMG (UNIQUE)
GL_ZE_vs_FMG_B2B1_enr
GL_ZE_vs_FMG_B3B2_enr
GL_ZE_vs_FMG_B4B3_enr
GL_ZE_vs_FMG_B5B4_enr
GL_ZE_vs_FMG_B6B5_enr
GL_ZE_vs_FMG_B7B6_enr
GL_ZE_vs_FMG_B8B7_enr
GL_ZE_vs_FMG_B9B8_enr
GL_ZE_vs_FMG_B10B9_enr

#PLOTS ZE vs FMG
{
  #ZE_vs_FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_vs_FMG_B2B1_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_vs_FMG_B2B1_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_vs_FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_vs_FMG_B3B2_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_vs_FMG_B3B2_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_vs_FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_vs_FMG_B5B4_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_vs_FMG_B5B4_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_vs_FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_vs_FMG_B6B5_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_vs_FMG_B6B5_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_vs_FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_vs_FMG_B7B6_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_vs_FMG_B7B6_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_vs_FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_vs_FMG_B8B7_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_vs_FMG_B8B7_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_vs_FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_vs_FMG_B9B8_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_vs_FMG_B9B8_enr, enrichment = "mapman", clusterColor = "#1B75BC")
  
  #ZE_vs_FMG_B2 GO and Mapman Treemap
  plotEnrichedTreemap(GL_ZE_vs_FMG_B10B9_enr, enrichment = "go", namespace = "none")
  plotEnrichedTreemap(GL_ZE_vs_FMG_B10B9_enr, enrichment = "mapman", clusterColor = "#1B75BC")
}





