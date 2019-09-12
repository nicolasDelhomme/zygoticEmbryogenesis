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
