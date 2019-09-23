library(tidyverse)

## Table of only "Pab" identities (FMG and ZE included)

### sed -i 's/_MIR1//g' sample_info.csv
### sed -i 's/_MIR2//g' sample_info.csv
### sed -i 's/_TR1//g' U.Egertsdotter_18_01_sample_info.csv
### sed -i 's/_TR2//g' U.Egertsdotter_18_01_sample_info.csv

### separate .csv files that are selecting out only the relevant column of "Pab" identities
pooled1 <- select(sample_info, 3)
pooled2 <- select(U.Egertsdotter_18_01_sample_info, 2)

### Unifying the table, removing duplicates and keeping differences.
unification <- union(pooled1, pooled2)
unification

write_csv(unification, path = "doc/testmerge.csv")


### if we want to remove FMG too, use: 
### sed -i 's/*_FMG//g' sample_info.csv 
### and 
### sed -i 's/*_FMG//g' U.Egertsdotter_18_01_sample_info.csv






### if we want to keep other data, remove first columns

pooled1 <- select(sample_info, 3,4,5,6,7)

translocate <- select(U.Egertsdotter_18_01_sample_info, 1)
pooled2 <- select(U.Egertsdotter_18_01_sample_info, 2,3,4)
pooled2 <- bind_cols(pooled2, translocate)
pooled1
pooled2





unification <- full_join(pooled1, pooled2)
unification
##unification <- inner_join(select(sample_info, 3,4,5,6,7), select(U.Egertsdotter_18_01_sample_info, 2,3,4))



write_csv(unification, path = "doc/testmerge_complete_v1.csv")



## sample PabB3C6_FMG and PabB3C6 are missing info
## sample PabB3C6_FMG time B3, sample id B3C6, tissue FMG, replicate B4-FMG - need to ammend this.

trashingInvalid <- testmerge_complete
print(trashingInvalid)

test <- slice(trashingInvalid,1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,
              32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58)
print(test)
write_csv(test, path = "doc/testmerge_complete_v2.csv")


arrangedV2 <- select(testmerge_complete_v2, 1,2,3,4,5,6,7)
arrangingNGI <- select(testmerge_complete_v2, 8)
print(arrangingNGI)
print(arrangedV2)
testjoin <- bind_cols(arrangingNGI, arrangedV2)
print(testjoin)
write_csv(testjoin, path = "doc/testmerge_complete_v3.csv")


#####adding Somatic and 29Seed information


filelistSomatic <- list.files("/mnt/picea/projects/spruce/uegertsdotter/SE-germinants/salmon", 
                                 recursive = TRUE, 
                                 pattern = "quant.sf",
                                 full.names = TRUE)
print(filelistSomatic)
filelist29Seed <- list.files("/mnt/picea/projects/spruce/uegertsdotter/29_Spruce_Seeds_Project/salmon", 
                                recursive = TRUE, 
                                pattern = "quant.sf",
                                full.names = TRUE)
print(filelist29Seed)
print(seeds.samples)

sampleseeds <- seeds.samples %>%
  rename(NGI.ID = SciLifeID,
         User.ID = SampleName)


unification <- testmerge_complete_v3

unification2 <- full_join(unification, sampleseeds)
print(unification2)

sampleSEGerm <- somaticEmbryogenesisGerminants.samples %>%
  rename(NGI.ID = SciLifeID,
         User.ID = SampleName)

unification3 <- full_join(unification2, sampleSEGerm)
print(unification3)


sampleSE <- somaticEmbryogenesis.samples %>%
  rename(NGI.ID = ScilifeID,
         User.ID = SubmittedID,
         Tissue = Description)
##replace tissue name with simple name, depending on the stage number. Then, remove stages column.
sampleSE

sampleSEWorkingFile <- str_replace(sampleSE$Tissue,".*Non-Embryogenic","Non-EMB") %>%
  str_replace(".*Proliferation.*","PEM") %>%
  str_replace(".*Early stage maturing embryos.*","DKM") %>%
  str_replace(".*Start maturation.*","SM") %>%
  str_replace(".*Late stage maturing.*","LSM") %>%
  str_replace(".*Start dessication.*","SD") %>%
  str_replace(".*End of dessication.*","ED") %>%
  str_replace(".*oot.*","RO") %>%
  str_replace(".*lanting.*","PLS") %>%
  str_replace(".*Proliferation.*","PEM") %>%
  str_replace(".*Start-dess.*","SD") %>%
  str_replace(".*Start-Mat.*","DKM") %>%
  str_replace(".*Start-mat.*","LSM") %>%
  str_replace(".*Mid-stage-mat.*","SM") %>%
  str_replace(".*Dessication.*","ED")
  
sampleSEWorkingFile
sampleSE["Tissue"] <- sampleSEWorkingFile
sampleSE <- subset(sampleSE, select = -c(Stages))
sampleSE
sampleSEWorkingFile <- str_replace(sampleSE$NGI.ID,".*XX_","") %>%
  str_replace("_d.*","")
sampleSE["NGI.ID"] <- sampleSEWorkingFile
str_which(sampleSE$NGI.ID,"_L002")
sampleSE <- sampleSE[-c(str_which(sampleSE$NGI.ID,"_L002")),]
sampleSEWorkingFile <- str_replace(sampleSE$NGI.ID,"_S.*","")
sampleSE["NGI.ID"] <- sampleSEWorkingFile

sampleSE
##split P464_2 from P464_3
##split the ones that have B
##remerge with _2
##remerge with _3

sampleSE.2 <- sampleSE[grep("_2", sampleSE$NGI.ID), ]
sampleSE.2b <- sampleSE.2[grep("B", sampleSE.2$NGI.ID), ]
sampleSE.2b
sampleSE.2brestruct <- sampleSE.2b[grep("7B", sampleSE.2b$NGI.ID), ]
sampleSE.2b <- sampleSE.2b[grep("7B", sampleSE.2b$NGI.ID, invert = TRUE), ]



sampleSE.2b <- full_join(sampleSE.2b, sampleSE.2brestruct)
sampleSE.2b



sampleSE.2 <- sampleSE.2[grep("B", sampleSE.2$NGI.ID, invert = TRUE), ]
sampleSE.2
sampleSE.3 <- sampleSE[grep("_3", sampleSE$NGI.ID), ]
sampleSE.3
##order is _2, _2 and B, then _3

sampleSEmerge <- full_join(sampleSE.2, sampleSE.2b)
sampleSEmerge <- full_join(sampleSEmerge, sampleSE.3)
sampleSE <- sampleSEmerge

unification4 <- full_join(unification3, sampleSE)
print(unification4)


time <- unification4$Time
time
tissue <- unification4$Tissue
tissue
ngi <- unification4$NGI.ID
experiment <- ngi
experiment <- str_replace(ngi,"P11562.*","Zygotic")
experiment <- str_replace(experiment, "P1790.*","29Seed")
experiment <- str_replace(experiment, "P2228.*","Somatic Embryogenesis Germinants")
experiment <- str_replace(experiment, "P464.*","Somatic Embryogenesis")
experiment <- str_replace(experiment, "P7614.*","Somatic Embryogenesis")
experiment

unification4["Experiment"] <- experiment
unification4

write_csv(unification4, path = "doc/testmerge_complete_v5.csv")
















##old code
NGI.ID <- str_replace(filelistSomatic,".*XX_","")
NGI.ID <- str_replace(NGI.ID,"_sort.*","")
NGI.ID <- unique(NGI.ID)
print(NGI.ID)

Tissue <- (rep(c("S"),times=19))
Time <- (rep(c("U"),times=19))
User.ID <- (rep(c("NA"),times=19))
Sample.ID <- (rep(c("NA"),times=19))
Replicate <- (rep(c("NA"),times=19))

somaticAdd <- data.frame(NGI.ID, User.ID, Time, Sample.ID, Tissue, Replicate)

NGI.ID <- str_replace(filelist29Seed,".*XX_","")
NGI.ID <- str_replace(NGI.ID,"_sort.*","")
NGI.ID <- unique(NGI.ID)
print(NGI.ID)

Tissue <- (rep(c("29Seed"),times=12))
Time <- (rep(c("U"),times=12))
User.ID <- (rep(c("NA"),times=12))
Sample.ID <- (rep(c("NA"),times=12))
Replicate <- (rep(c("NA"),times=12))

seedAdd <- data.frame(NGI.ID, User.ID, Time, Sample.ID, Tissue, Replicate)


somaticAdd
seedAdd



unification2 <- bind_rows(testmerge_complete_v3,somaticAdd)
unification2

unification3 <- bind_rows(unification2,seedAdd)
unification3

write_csv(unification3, path = "doc/testmerge_complete_v4.csv")

