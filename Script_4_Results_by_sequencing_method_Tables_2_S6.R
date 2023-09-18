### ========================== SCRIPT 4 - COMPARISON OF SEQUENCING METHODS ========================== ###

### SCRIPT:          SCRIPT 4 - COMPARISON OF SEQUENCING METHODS
##
## DESCRIPTION:      Script to generate/run
##                    - Table 2
##                    - Table S6
##
## AUTHORS:          Johanne Spreckels
##
## NOTES:            To use this script, set correct paths to the following folders/files:
##                    - MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_DATA/
##                    - MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/
##
## DATE OF CREATION: Script generated for publication in September 2023


### CONTENTS OF THIS FILE
## 0. IMPORT DATA
## 1. SELECT SAMPLES OF INTEREST
## 2. SUBSET SEQUENCING DATA BY TAXONOMIC LEVEL
## 3. CALCULATE GENERA/SPECIES/ASV RICHNESS
## 4. TABLE 2A: SUMMARY STATISTICS FOR MOCK SAMPLES
## 5. TABLE 2B: SUMMARY STATISTICS FOR NEGATIVE CONTROLS
## 6. TABLE 2C: SUMMARY STATISTICS FOR MILK SAMPLES
## 7. TABLE S6: ASVs DETECTED IN NEGATIVE CONTROLS


### ===== 0. IMPORT DATA ===== ###

## set working directory
setwd("MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_DATA/")

## import Data_sequencing_method_comparison_n136.txt
data <- read.table("Data_sequencing_method_comparison_n136.txt", header=T, sep="\t")


### ===== 1. SELECT SAMPLES OF INTEREST ===== ###

## set alias as rownames
rownames(data) <- as.factor(as.character(paste0(data$sequencing_method, "_", data$subjectId)))
rownames(data)[grep("Mock", rownames(data))] <- c("16S-Mock", "16S-ITS-23S-Mock")

## select by sample type
mock    <- data[data$sample_set=="method_test_sample" & data$sample_type=="Mock",]             #2 samples
negctrl <- data[data$sample_set=="method_test_sample" & data$sample_type=="Negative_control",] #7 samples
milk    <- data[data$sample_set=="method_test_sample" & data$sample_type=="Milk",]             #15 samples


### ===== 2. SUBSET SEQUENCING DATA BY TAXONOMIC LEVEL ===== ###

## function to subset data frame by taxonomic level (genus levels)
subset.by.taxlevel <- function(inputdata, pattern){ 
  colnames(inputdata) <- gsub(pattern, "", colnames(inputdata)) #shorten colnames after selected taxonomic level
  taxon <- t(rowsum(t(inputdata), group = colnames(inputdata), na.rm = T)) #combine reads from columns with the same name in one column
  taxon <- taxon[,colSums(taxon)>0] #remove columns where colSums==0
  return(taxon)
}

## subset mock data by taxonomic level
mock_genus <- subset.by.taxlevel(inputdata = mock[,22:ncol(mock)], pattern = "_s__.*") #2x8 genera

## subset negctrl data by taxonomic level
negctrl_genus <- subset.by.taxlevel(inputdata = negctrl[,22:ncol(negctrl)], pattern = "_s__.*") #7x21 genera

## subset milk data by taxonomic level
milk_genus <- subset.by.taxlevel(inputdata = milk[,22:ncol(milk)], pattern = "_s__.*") #15x65 genera


## create data frames for ASV level

## mock
mock_bact <- mock[,22:ncol(mock)]
mock_ASV <- mock_bact[,colSums(mock_bact)>0] #2x23 ASVs

## negctrl
negctrl_bact <- negctrl[,22:ncol(mock)]
negctrl_ASV <- negctrl_bact[,colSums(negctrl_bact)>0] #7x25 ASVs

## milk
milk_bact <- milk[,22:ncol(mock)]
milk_ASV <- milk_bact[,colSums(milk_bact)>0] #15x270 ASVs


### ===== 3. CALCULATE GENERA/SPECIES/ASV RICHNESS ===== ###

library(vegan)

#for each sample type, calculate richness for genus level and for ASV level - only taking into account bacteria that could be classified on the given level

### MOCK

## calculate richness
mock_genus_richness <- as.data.frame(specnumber(mock_genus[,grep("g__NA", colnames(mock_genus), invert=T)]))
mock_ASV_richness <- as.data.frame(specnumber(mock_ASV[,grep("ASV__NA", colnames(mock_ASV), invert=T)]))

## combine in 1 data frame
mock_richness <- as.data.frame(cbind(mock_genus_richness, mock_ASV_richness))
mock_richness$sample_name <- rownames(mock_richness)
mock_richness$sequencing_method <- c("16S_V3V4", "16S-ITS-23S")
colnames(mock_richness)[1:2] <- c("genera_richness", "ASV_richness")
mock_richness <- mock_richness[,c(3:4,1:2)]
mock_richness


### NEGATIVE CONTROLS

## calculate richness
negctrl_genus_richness <- as.data.frame(specnumber(negctrl_genus[,grep("g__NA", colnames(negctrl_genus), invert=T)]))
negctrl_ASV_richness <- as.data.frame(specnumber(negctrl_ASV[,grep("ASV__NA", colnames(negctrl_ASV), invert=T)]))

## combine in 1 data frame
negctrl_richness <- as.data.frame(cbind(negctrl_genus_richness, negctrl_ASV_richness))
negctrl_richness$sample_name <- rownames(negctrl_richness)
negctrl_richness$sequencing_method <- c(rep("16S_V3V4",4), rep("16S-ITS-23S",3))
colnames(negctrl_richness)[1:2] <- c("genera_richness", "ASV_richness")
negctrl_richness <- negctrl_richness[,c(3:4,1:2)]
negctrl_richness[negctrl_richness$sample_name=="16S_V3V4_Isol-NC-1",3:4] <- NA #the negative control 16S_V3V4_Isol-NC-1 had 0 reads, so I set the richness to NA instead of 0
negctrl_richness


### MILK

## calculate richness
milk_genus_richness <- as.data.frame(specnumber(milk_genus[,grep("g__NA", colnames(milk_genus), invert=T)]))
milk_ASV_richness <- as.data.frame(specnumber(milk_ASV[,grep("ASV__NA", colnames(milk_ASV), invert=T)]))

## combine in 1 data frame
milk_richness <- as.data.frame(cbind(milk_genus_richness, milk_ASV_richness))
milk_richness$sample_name <- rownames(milk_richness)
milk_richness$sequencing_method <- c(rep("16S_V3V4",9), rep("16S-ITS-23S",6))
milk_richness <- milk_richness[order(milk_richness$sample_name),]
colnames(milk_richness)[1:2] <- c("genera_richness", "ASV_richness")
milk_richness <- milk_richness[,c(3:4,1:2)]
milk_richness


# add richness info to original data frames
mock_meta <- merge(mock[,c(1:21)], mock_richness[,-c(1:2)], by="row.names")
mock_meta <- mock_meta[,-1]

negctrl_meta <- merge(negctrl[,c(1:21)], negctrl_richness[,-c(1:2)], by="row.names")
negctrl_meta <- negctrl_meta[,-1]

milk_meta <- merge(milk[,c(1:21)], milk_richness[,-c(1:2)], by="row.names")
milk_meta <- milk_meta[,-1]


### ===== 4. TABLE 2A: SUMMARY STATISTICS FOR MOCK SAMPLES ===== ###

# split mock_meta data frame by sequencing method
mock_meta_16S23S <- mock_meta[mock_meta$sequencing_method=="16S-ITS-23S",]
mock_meta_16S_V3V4 <- mock_meta[mock_meta$sequencing_method=="16S_V3V4",]

# function to retrieve values
getSumStats_mock <- function(inputdata){
  SumStats_mock <- c(nrow(inputdata),
                      inputdata$n_reads_total,
                      inputdata$ASV_richness,
                      inputdata$genera_richness)
  return(SumStats_mock)
}

SumStats_mock <- as.data.frame(cbind(getSumStats_mock(mock_meta_16S23S),
                                     getSumStats_mock(mock_meta_16S_V3V4)))
SumStats_mock[5,] <- c(rep("Bacterial mock community",2))
SumStats_mock[6,] <- c("16S-ITS-23S","16S_V3V4")
SumStats_mock$info <- c("n (samples)", "n (reads)", "n (ASVs)", "n (classified genera)", "Sample", "Sequencing method")
SumStats_mock <- SumStats_mock[c(5:6,1:4), c(3,1,2)]
SumStats_mock
# write.table(SumStats_mock, file="MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/Table_2A.txt", sep="\t", row.names=F, quote=F)


### ===== 5. TABLE 2B: SUMMARY STATISTICS FOR NEGATIVE CONTROLS ===== ###

# split negctrl_meta data frame by sequencing method and type of negative control (isollation vs library prep negative control)
isol_negctrl_meta <- negctrl_meta[grep("Isol-NC", negctrl_meta$sample_ID),]
isol_negctrl_meta_16S23S <- isol_negctrl_meta[isol_negctrl_meta$sequencing_method=="16S-ITS-23S",]
isol_negctrl_meta_16S_V3V4 <- isol_negctrl_meta[isol_negctrl_meta$sequencing_method=="16S_V3V4",]

libprep_negctrl_meta <- negctrl_meta[grep("LibPrep-NC", negctrl_meta$sample_ID),]
libprep_negctrl_meta_16S23S <- libprep_negctrl_meta[libprep_negctrl_meta$sequencing_method=="16S-ITS-23S",]
libprep_negctrl_meta_16S_V3V4 <- libprep_negctrl_meta[libprep_negctrl_meta$sequencing_method=="16S_V3V4",]

# function to retrieve median (range)
getSumStats_negctrl <- function(inputdata){
  SumStats_negctrl <- c(nrow(inputdata),
                      paste0(median(inputdata$n_reads_total), " (", min(inputdata$n_reads_total), " - ", max(inputdata$n_reads_total), ")"),
                      paste0(median(inputdata$ASV_richness, na.rm=T), " (", min(inputdata$ASV_richness, na.rm=T), " - ", max(inputdata$ASV_richness, na.rm=T), ")"), #note: remove the negctrl with 0 reads here
                      paste0(median(inputdata$genera_richness, na.rm=T), " (", min(inputdata$genera_richness, na.rm=T), " - ", max(inputdata$genera_richness, na.rm=T), ")")) #note: remove the negctrl with 0 reads here
  return(SumStats_negctrl)
}

SumStats_negctrl <- as.data.frame(cbind(getSumStats_negctrl(isol_negctrl_meta_16S23S),
                                        getSumStats_negctrl(libprep_negctrl_meta_16S23S),
                                        getSumStats_negctrl(isol_negctrl_meta_16S_V3V4),
                                        getSumStats_negctrl(libprep_negctrl_meta_16S_V3V4)))
for (i in 1:ncol(SumStats_negctrl)){SumStats_negctrl[,i] <- as.character(as.factor(SumStats_negctrl[,i]))}
SumStats_negctrl[5,] <- c("DNA isolation negative controls", "Library preparation negative control", "DNA isolation negative controls", "Library preparation negative control")
SumStats_negctrl[6,] <- c(rep("16S-ITS-23S",2), rep("16S_V3V4",2))
SumStats_negctrl$info <- c("n (samples)", "n (reads)", "n (ASVs)", "n (classified genera)", "Samples", "Sequencing method")
SumStats_negctrl <- SumStats_negctrl[c(5:6,1:4), c(5,1:4)]
SumStats_negctrl
# write.table(SumStats_negctrl, file="MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/Table_2B.txt", sep="\t", row.names=F, quote=F)


### ===== 6. TABLE 2C: SUMMARY STATISTICS FOR MILK SAMPLES ===== ###

## split milk_meta data frame by DNA isolation method
milk_meta_16S23S <- milk_meta[milk_meta$sequencing_method=="16S-ITS-23S",]
milk_meta_16S_V3V4 <- milk_meta[milk_meta$sequencing_method=="16S_V3V4",]

## function to retrieve median (range)
getSumStats_milk <- function(inputdata){
  SumStats_milk <- c(nrow(inputdata),
                     paste0(median(inputdata$n_reads_total), " (", min(inputdata$n_reads_total), " - ", max(inputdata$n_reads_total), ")"),
                     paste0(median(inputdata$ASV_richness), " (", min(inputdata$ASV_richness), " - ", max(inputdata$ASV_richness), ")"),
                     paste0(median(inputdata$genera_richness), " (", min(inputdata$genera_richness), " - ", max(inputdata$genera_richness), ")"))
  
  names(SumStats_milk) <- c("n (samples)", "n (reads)", "n (ASVs)", "n (classified genera)")
  
  return(SumStats_milk)
}

SumStats_milk <- as.data.frame(cbind(getSumStats_milk(milk_meta_16S23S),
                                     getSumStats_milk(milk_meta_16S_V3V4)))
SumStats_milk$info <- rownames(SumStats_milk)

# function to get p and FDR
getStats_milk <- function(inputdata, columns_of_interest){
  p <- c()
  
  for (i in columns_of_interest){
    print(colnames(inputdata)[i])
    wt <- wilcox.test(inputdata[inputdata$sequencing_method=="16S-ITS-23S",i], inputdata[inputdata$sequencing_method=="16S_V3V4",i], alternative="two.sided", paired=F, exact=T)
    p <- c(p, wt$p.value)
  }
  
  FDR <- p.adjust(p, method = "BH")
  results <- data.frame(p=round(p,4), FDR=round(FDR,3))
  rownames(results)<- c("n (reads)", "n (ASVs)", "n (classified genera)")
  
  return(results)
}

Stats_milk <- getStats_milk(milk_meta, columns_of_interest = c(21,23,22))
Stats_milk$info <- rownames(Stats_milk)

## add statistics to milk summary stats table
SumStats_milk2 <- left_join(SumStats_milk[,c(3,1:2)], Stats_milk, by="info")

## clean up final milk summary stats table
for (i in 1:3){SumStats_milk2[,i] <- as.character(as.factor(SumStats_milk2[,i]))}

SumStats_milk2[5,] <- c("Samples", "3 milk samples in duplicate", "3 milk samples in triplicate", NA, NA)
SumStats_milk2[6,] <- c("Sequencing method", "16S-ITS-23S", "16S_V3V4", NA, NA)

SumStats_milk2 <- SumStats_milk2[c(6,5,1:4),]
SumStats_milk2
# write.table(SumStats_milk2, file="MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/Table_2C.txt", sep="\t", row.names=F, quote=F)


### ===== 7. TABLE S6: ASVs DETECTED IN NEGATIVE CONTROLS ===== ###

library(dplyr)

## split negctrl_ASV dataframe by sequencing method and type of negative control (isolation vs libprep negative control)
#16S-23S
Seq_16S23S_negctrl_ASV <- negctrl_ASV[grep("16S-ITS-23S_Isol-NC", rownames(negctrl_ASV)),]
Seq_16S23S_libprep_negctrl_ASV <- negctrl_ASV[grep("16S-ITS-23S_LibPrep-NC", rownames(negctrl_ASV)),]

#16S_V3V4
Seq_16S_V3V4_negctrl_ASV <- negctrl_ASV[grep("16S_V3V4_Isol-NC", rownames(negctrl_ASV)),]
Seq_16S_V3V4_libprep_negctrl_ASV <- negctrl_ASV[grep("16S_V3V4_LibPrep-NC", rownames(negctrl_ASV)),]

## remove columns with 0 reads for the selected negative controls
#16S-ITS-23S
Seq_16S23S_negctrl_ASV <- Seq_16S23S_negctrl_ASV[,colSums(Seq_16S23S_negctrl_ASV)!=0] #2 different ASVs in 16S-23S isolation controls
Seq_16S23S_libprep_negctrl_ASV_tmp <- as.data.frame(Seq_16S23S_libprep_negctrl_ASV[,colSums(Seq_16S23S_libprep_negctrl_ASV)!=0]) #1 ASV in 16S-23S libprep control
rownames(Seq_16S23S_libprep_negctrl_ASV_tmp) <- rownames(Seq_16S23S_libprep_negctrl_ASV)
colnames(Seq_16S23S_libprep_negctrl_ASV_tmp) <- colnames(Seq_16S23S_libprep_negctrl_ASV)[24]
Seq_16S23S_libprep_negctrl_ASV <- Seq_16S23S_libprep_negctrl_ASV_tmp

#16S_V3V4
Seq_16S_V3V4_negctrl_ASV <- Seq_16S_V3V4_negctrl_ASV[,colSums(Seq_16S_V3V4_negctrl_ASV)!=0] #19 different ASVs in 16S isolation controls
Seq_16S_V3V4_libprep_negctrl_ASV <- Seq_16S_V3V4_libprep_negctrl_ASV[,colSums(Seq_16S_V3V4_libprep_negctrl_ASV)!=0] #4 different ASVs in 16S libprep control


## save supplementary table showing which bacteria / ASVs are detected in negative controls (1 table for all sequencing methods)

# ensure rownames and sample_name columns match before merging
df_negctrl_ASV <- negctrl_ASV
df_negctrl_ASV$sample_name <- rownames(df_negctrl_ASV)
rownames(negctrl_meta) <- as.factor(as.character(paste0(negctrl_meta$sequencing_method, "_", negctrl_meta$subjectId)))
negctrl_meta$sample_name <- rownames(negctrl_meta)

# create joined table
df_negctrl_ASV2 <- left_join(df_negctrl_ASV, negctrl_meta[,c(24,21)], by="sample_name") #add total number of reads to the table
df_negctrl_ASV2 <- df_negctrl_ASV2[,c(26,27,1:25)] #show sample name and read count columns first
df_negctrl_ASV2[,3:ncol(df_negctrl_ASV2)] <- df_negctrl_ASV2[,3:ncol(df_negctrl_ASV2)] / df_negctrl_ASV2$n_reads_total #calculate relative abundances
df_negctrl_ASV2[df_negctrl_ASV2$n_reads_total==0,3:27] <- NA
colnames(df_negctrl_ASV2)[1] <- "sample"
# write.table(df_negctrl_ASV2, file="MY_PATH_TO_SEQUENCING_METHOD_COMPARISON_RESULTS/Table_S6.txt", sep="\t", row.names=F, quote=F)












