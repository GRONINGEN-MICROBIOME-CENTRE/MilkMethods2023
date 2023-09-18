### ========================== SCRIPT 1 - COMPARISON OF DNA ISOLATION METHODS ========================== ###

### SCRIPT:          SCRIPT 1 - RESULTS BY DNA ISOLATION METHOD
##
## DESCRIPTION:      Script to generate/run
##                    - Table 1
##                    - Table S2
##                    - Table S4
##                    - the Wilcoxon test to compare read counts between negative controls and milk samples
##                    - the correlation between read count or DNA yield and the percentage of reads shared with negative controls
##
## AUTHORS:          Johanne Spreckels
##
## NOTES:            To use this script, set correct paths to the following folders/files:
##                    - MY_PATH_TO_DNA_ISOLATION_COMPARISON_DATA/Data_DNA_isolation_comparison_n42.txt
##                    - MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/
##
## DATE OF CREATION: Script generated for publication in September 2023

### CONTENTS OF THIS FILE
## 0. IMPORT DATA
## 1. SUBSET SEQUENCING DATA BY SAMPLE TYPE
## 2. SUBSET SEQUENCING DATA BY TAXONOMIC LEVEL
## 3. CALCULATE GENERA AND ASV RICHNESS
## 4. DETECTING POTENTIAL CONTAMINATION: ASVs IN MILK SAMPLES SHARED WITH NEGATIVE CONTROLS AND UNIQUE TO MILK (TABLE S2, TABLE S4)
## 5. MILK SAMPLES ONLY: CALCULATE RICHNESS AFTER EXCLUDING CONTAMINATION
## 6. MOCK AND MILK SAMPLES ONLY: CALCULATE PERCENTAGE OF READS CLASSIFIED ON SPECIES LEVEL
## 7. TABLE 1A: SUMMARY STATISTICS FOR MOCK SAMPLES
## 8. TABLE 1B: SUMMARY STATISTICS FOR NEGATIVE CONTROLS
## 9. TABLE 1C: SUMMARY STATISTICS FOR MILK SAMPLES
## 10. COMPARISON OF READ COUNT BETWEEN MILK AND NEGATIVE CONTROLS
## 11. CORRELATION OF READ COUNT AND DNA YIELD WITH PERCENTAGE OF CONTAMINATION IN MILK SAMPLES


### ===== 0. IMPORT DATA ===== ###

## import Data_DNA_isolation_comparison_n42.txt data file
# note that this data file contains ASVs, which are not assigned on certain taxonomic levels, incl. on kingdom and phylum level
data <- read.table("MY_PATH_TO_DNA_ISOLATION_COMPARISON_DATA/Data_DNA_isolation_comparison_n42.txt", header=T, sep="\t", stringsAsFactors=T)


### ===== 1. SUBSET SEQUENCING DATA BY SAMPLE TYPE ===== ###

## set sample_ids as rownames
rownames(data) <- data$sample_ID

## subset by sample type
mock <- data[data$sample_type=="Mock",] #5 samples
milk <- data[data$sample_type=="Milk",] #27 samples
negctrl <- data[data$sample_type=="Negative_control",] #10 samples; note that 1 PSK-isolated negative control (BM_N1_B_NC_1) has 0 reads!


### ===== 2. SUBSET SEQUENCING DATA BY TAXONOMIC LEVEL ===== ###

## function to subset data frame by taxonomic level
subset.by.taxlevel <- function(inputdata, pattern){ 
  colnames(inputdata) <- gsub(pattern, "", colnames(inputdata)) #shorten colnames after selected taxonomic level
  taxon <- t(rowsum(t(inputdata), group = colnames(inputdata), na.rm = T)) #combine reads from columns with the same name in one column
  taxon <- taxon[,colSums(taxon)>0] #remove columns where colSums==0
  return(taxon)
}

## subset mock data by taxonomic level
mock_genus <- as.data.frame(subset.by.taxlevel(inputdata = mock[,19:ncol(mock)], pattern = "_s__.*")) #5x9 genera
mock_species <- as.data.frame(subset.by.taxlevel(inputdata = mock[,19:ncol(mock)], pattern = "_ASV__.*")) #5x11 species

## subset negctrl data by taxonomic level
negctrl_genus <- as.data.frame(subset.by.taxlevel(inputdata = negctrl[,19:ncol(negctrl)], pattern = "_s__.*")) #10x75 genera

## subset milk data by taxonomic level
milk_genus <- as.data.frame(subset.by.taxlevel(inputdata = milk[,19:ncol(milk)], pattern = "_s__.*")) #27x167 genera


## create data frames for ASV level

## mock
mock_bact <- mock[,19:ncol(mock)]
mock_ASV <- mock_bact[,colSums(mock_bact)>0] #5x17 ASVs

## negctrl
negctrl_bact <- negctrl[,19:ncol(negctrl)]
negctrl_ASV <- negctrl_bact[,colSums(negctrl_bact)>0] #10x101 ASVs

## milk
milk_bact <- milk[,19:ncol(milk)]
milk_ASV <- milk_bact[,colSums(milk_bact)>0] #27x412 ASVs


### ===== 3. CALCULATE GENERA AND ASV RICHNESS ===== ###

library(vegan)

#for each sample type, calculate richness for genus level and for ASV level - only taking into account bacteria that could be classified on the given level

### MOCK

## calculate richness
mock_genus_richness <- as.data.frame(specnumber(mock_genus[,grep("g__NA", colnames(mock_genus), invert=T)]))
mock_ASV_richness <- as.data.frame(specnumber(mock_ASV[,grep("ASV__NA", colnames(mock_ASV), invert=T)]))

## combine in 1 data frame
mock_richness <- as.data.frame(cbind(mock_genus_richness, mock_ASV_richness))
mock_richness$sample_ID <- rownames(mock_richness)
mock_richness <- mock_richness[order(mock_richness$sample_ID),]
mock_richness$DNA_isolation_method <- c("MilkDNAKit", "PowerSoilProKit", "MagMAXKit", "FastStoolMiniKit", "DNA_Mock")
colnames(mock_richness)[1:2] <- c("genera_richness", "ASV_richness")
mock_richness <- mock_richness[,c(3:4,1:2)]
mock_richness


### NEGATIVE CONTROLS

## calculate richness
negctrl_genus_richness <- as.data.frame(specnumber(negctrl_genus[,grep("g__NA", colnames(negctrl_genus), invert=T)]))
negctrl_ASV_richness <- as.data.frame(specnumber(negctrl_ASV[,grep("ASV__NA", colnames(negctrl_ASV), invert=T)]))

## combine in 1 data frame
negctrl_richness <- as.data.frame(cbind(negctrl_genus_richness, negctrl_ASV_richness))
negctrl_richness$sample_ID <- rownames(negctrl_richness)
negctrl_richness <- negctrl_richness[order(negctrl_richness$sample_ID),]
negctrl_richness$DNA_isolation_method <- c(rep("MilkDNAKit",2), rep("PowerSoilProKit",3), rep("MagMAXKit",2), rep("FastStoolMiniKit",2), "Library_prep_negctrl")
colnames(negctrl_richness)[1:2] <- c("genera_richness", "ASV_richness")
negctrl_richness <- negctrl_richness[,c(3:4,1:2)]
negctrl_richness[negctrl_richness$sample_ID=="BM_N1_B_NC_1",3:4] <- NA #the negative control BM_N1_B_NC_1 had 0 reads, so I set the richness to NA instead of 0
negctrl_richness


### MILK

## calculate richness
milk_genus_richness <- as.data.frame(specnumber(milk_genus[,grep("g__NA", colnames(milk_genus), invert=T)]))
milk_ASV_richness <- as.data.frame(specnumber(milk_ASV[,grep("ASV__NA", colnames(milk_ASV), invert=T)]))

## combine in 1 data frame
milk_richness <- as.data.frame(cbind(milk_genus_richness, milk_ASV_richness))
milk_richness$sample_ID <- rownames(milk_richness)
milk_richness <- milk_richness[order(milk_richness$sample_ID),]
milk_richness$DNA_isolation_method <- c(rep("MilkDNAKit",6), rep("PowerSoilProKit",9), rep("MagMAXKit",6), rep("FastStoolMiniKit",6))
colnames(milk_richness)[1:2] <- c("genera_richness", "ASV_richness")
milk_richness <- milk_richness[,c(3:4,1:2)]
milk_richness


# add richness info to original data frames
mock_meta <- merge(mock[,1:18], mock_richness[,-c(1:2)], by="row.names")
mock_meta <- mock_meta[,-1]

negctrl_meta <- merge(negctrl[,1:18], negctrl_richness[,-c(1:2)], by="row.names")
negctrl_meta <- negctrl_meta[,-1]

milk_meta <- merge(milk[,1:18], milk_richness[,-c(1:2)], by="row.names")
milk_meta <- milk_meta[,-1]


### ===== 4. DETECTING POTENTIAL CONTAMINATION: ASVs IN MILK SAMPLES SHARED WITH NEGATIVE CONTROLS AND UNIQUE TO MILK (TABLE S2, TABLE S4) ===== ###

library(dplyr)

### NEGATIVE CONTROLS -> identify ASVs in negative controls

## split negctrl_ASV dataframe by DNA isolation method
MDK_negctrl_ASV <- negctrl_ASV[grep("_A_", rownames(negctrl_ASV)),]
PSK_negctrl_ASV <- negctrl_ASV[grep("_B_", rownames(negctrl_ASV)),]
MXK_negctrl_ASV <- negctrl_ASV[grep("_D_", rownames(negctrl_ASV)),]
FSK_negctrl_ASV <- negctrl_ASV[grep("_E_", rownames(negctrl_ASV)),]
Seq_negctrl_ASV <- negctrl_ASV[grep("Library_prep_NTC", rownames(negctrl_ASV)),]

## remove columns with 0 reads for the selected negative controls
MDK_negctrl_ASV <- MDK_negctrl_ASV[,colSums(MDK_negctrl_ASV)!=0] #33 different ASVs in MDK isolation controls
PSK_negctrl_ASV <- PSK_negctrl_ASV[,colSums(PSK_negctrl_ASV)!=0] #19 different ASVs in MDK isolation controls
MXK_negctrl_ASV <- MXK_negctrl_ASV[,colSums(MXK_negctrl_ASV)!=0] # 9 different ASVs in MDK isolation controls
FSK_negctrl_ASV <- FSK_negctrl_ASV[,colSums(FSK_negctrl_ASV)!=0] #48 different ASVs in MDK isolation controls
Seq_negctrl_ASV <- Seq_negctrl_ASV[,colSums(Seq_negctrl_ASV)!=0] # 4 different ASVs in MDK isolation controls

## save supplementary table showing which bacteria / ASVs are detected in negative controls (1 table for all DNA isolation methods)
df_negctrl_ASV <- negctrl_ASV[order(rownames(negctrl_ASV)),]
df_negctrl_ASV$sample <- c("MilkDNAKit_NC1", "MilkDNAKit_NC2", "PowerSoilProKit_NC1", "PowerSoilProKit_NC2", "PowerSoilProKit_NC3", "MagMAXKit_NC1", "MagMAXKit_NC2", "FastDNAStoolMiniKit_NC1","FastDNAStoolMiniKit_NC2", "Library_preparation_NC1")
df_negctrl_ASV$sample_ID <- rownames(df_negctrl_ASV)
df_negctrl_ASV2 <- left_join(df_negctrl_ASV, negctrl_meta[,c(13,18)], by="sample_ID") #add total number of reads to the table
df_negctrl_ASV2 <- df_negctrl_ASV2[,c(102,104,1:101)] #show sample and read count columns first
df_negctrl_ASV2[,3:ncol(df_negctrl_ASV2)] <- df_negctrl_ASV2[,3:ncol(df_negctrl_ASV2)] / df_negctrl_ASV2$n_reads_dada2_final #calculate relative abundances
df_negctrl_ASV2[3,3:ncol(df_negctrl_ASV2)] <- NA #the negative control BM_N1_B_NC_1 had 0 reads, so I set the relative abundances to NA
df_negctrl_ASV2 <- df_negctrl_ASV2[c(8:9,1:7,10),] #resort rows as in other figures and tables
colnames(df_negctrl_ASV2)[2] <- "n_reads_total"
# write.table(df_negctrl_ASV2, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Table_S2.txt", sep="\t", row.names=F, quote=F)


### MILK SAMPLES -> check which milk ASVs are also detected in the corresponding isolation / library prep negative control(s)

## split milk_ASV dataframe by DNA isolation method
MDK_milk_ASV <- milk_ASV[grep("_A_", rownames(milk_ASV)),]
PSK_milk_ASV <- milk_ASV[grep("_B_", rownames(milk_ASV)),]
MXK_milk_ASV <- milk_ASV[grep("_D_", rownames(milk_ASV)),]
FSK_milk_ASV <- milk_ASV[grep("_E_", rownames(milk_ASV)),]

## remove columns with 0 reads for the selected milk samples
MDK_milk_ASV <- MDK_milk_ASV[,colSums(MDK_milk_ASV)!=0] #156 ASVs detected
PSK_milk_ASV <- PSK_milk_ASV[,colSums(PSK_milk_ASV)!=0] #138 ASVs detected
MXK_milk_ASV <- MXK_milk_ASV[,colSums(MXK_milk_ASV)!=0] #179 ASVs detected
FSK_milk_ASV <- FSK_milk_ASV[,colSums(FSK_milk_ASV)!=0] #191 ASVs detected

## identify ASVs shared between milk and their corresponding isolation negative controls
MDK_milk_ASV_isol_negctrl <- MDK_milk_ASV[,(colnames(MDK_milk_ASV) %in% colnames(MDK_negctrl_ASV))] #13 ASVs shared
PSK_milk_ASV_isol_negctrl <- PSK_milk_ASV[,(colnames(PSK_milk_ASV) %in% colnames(PSK_negctrl_ASV))] # 4 ASVs shared
MXK_milk_ASV_isol_negctrl <- MXK_milk_ASV[,(colnames(MXK_milk_ASV) %in% colnames(MXK_negctrl_ASV))] # 3 ASVs shared
FSK_milk_ASV_isol_negctrl <- FSK_milk_ASV[,(colnames(FSK_milk_ASV) %in% colnames(FSK_negctrl_ASV))] #17 ASVs shared

## identify ASVs shared between milk and the library preparation negative control
MDK_milk_ASV_Seq_negctrl <- MDK_milk_ASV[,(colnames(MDK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] #2 ASVs shared
PSK_milk_ASV_Seq_negctrl <- PSK_milk_ASV[,(colnames(PSK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] #0 ASVs shared
MXK_milk_ASV_Seq_negctrl <- MXK_milk_ASV[,(colnames(MXK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] #0 ASVs shared
FSK_milk_ASV_Seq_negctrl <- FSK_milk_ASV[,(colnames(FSK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] #2 ASVs shared

## identify ASVs shared between milk and any corresponding negative control (isolation and/or library prep negative control)
MDK_milk_ASV_isol_Seq_negctrl <- MDK_milk_ASV[,(colnames(MDK_milk_ASV) %in% colnames(MDK_negctrl_ASV)) | (colnames(MDK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] #15 ASVs shared
PSK_milk_ASV_isol_Seq_negctrl <- PSK_milk_ASV[,(colnames(PSK_milk_ASV) %in% colnames(PSK_negctrl_ASV)) | (colnames(PSK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] # 4 ASVs shared
MXK_milk_ASV_isol_Seq_negctrl <- MXK_milk_ASV[,(colnames(MXK_milk_ASV) %in% colnames(MXK_negctrl_ASV)) | (colnames(MXK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] # 3 ASVs shared
FSK_milk_ASV_isol_Seq_negctrl <- FSK_milk_ASV[,(colnames(FSK_milk_ASV) %in% colnames(FSK_negctrl_ASV)) | (colnames(FSK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] #18 ASVs shared

## calculate number of reads belonging to ASVs shared with corresponding negative controls
MDK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_negctrls <- rowSums(MDK_milk_ASV_isol_Seq_negctrl)
PSK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_negctrls <- rowSums(PSK_milk_ASV_isol_Seq_negctrl)
MXK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_negctrls <- rowSums(MXK_milk_ASV_isol_Seq_negctrl)
FSK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_negctrls <- rowSums(FSK_milk_ASV_isol_Seq_negctrl)

## calculate number of reads belonging to ASVs shared with corresponding isolation negative controls
MDK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_isolation_negctrls <- rowSums(MDK_milk_ASV_isol_negctrl)
PSK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_isolation_negctrls <- rowSums(PSK_milk_ASV_isol_negctrl)
MXK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_isolation_negctrls <- rowSums(MXK_milk_ASV_isol_negctrl)
FSK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_isolation_negctrls <- rowSums(FSK_milk_ASV_isol_negctrl)

## calculate number of reads belonging to ASVs shared with corresponding library preparation negative control
MDK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_library_preparation_negctrls <- rowSums(MDK_milk_ASV_Seq_negctrl)
PSK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_library_preparation_negctrls <- rowSums(PSK_milk_ASV_Seq_negctrl)
MXK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_library_preparation_negctrls <- rowSums(MXK_milk_ASV_Seq_negctrl)
FSK_milk_ASV_isol_Seq_negctrl$n_reads_shared_with_library_preparation_negctrls <- rowSums(FSK_milk_ASV_Seq_negctrl)

## add column showing sample_ID to tables
MDK_milk_ASV_isol_Seq_negctrl$sample_ID <- as.factor(as.character(rownames(MDK_milk_ASV_isol_Seq_negctrl)))
PSK_milk_ASV_isol_Seq_negctrl$sample_ID <- as.factor(as.character(rownames(PSK_milk_ASV_isol_Seq_negctrl)))
MXK_milk_ASV_isol_Seq_negctrl$sample_ID <- as.factor(as.character(rownames(MXK_milk_ASV_isol_Seq_negctrl)))
FSK_milk_ASV_isol_Seq_negctrl$sample_ID <- as.factor(as.character(rownames(FSK_milk_ASV_isol_Seq_negctrl)))

## add total read count to tables
MDK_milk_ASV_isol_Seq_negctrl_reads <- left_join(MDK_milk_ASV_isol_Seq_negctrl, milk[,c(13,18)], by="sample_ID")
PSK_milk_ASV_isol_Seq_negctrl_reads <- left_join(PSK_milk_ASV_isol_Seq_negctrl, milk[,c(13,18)], by="sample_ID")
MXK_milk_ASV_isol_Seq_negctrl_reads <- left_join(MXK_milk_ASV_isol_Seq_negctrl, milk[,c(13,18)], by="sample_ID")
FSK_milk_ASV_isol_Seq_negctrl_reads <- left_join(FSK_milk_ASV_isol_Seq_negctrl, milk[,c(13,18)], by="sample_ID")

## add percentage of reads belonging to ASVs shared with corresponding negative controls to table
MDK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_negctrls <- MDK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_negctrls / MDK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
PSK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_negctrls <- PSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_negctrls / PSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
MXK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_negctrls <- MXK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_negctrls / MXK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
FSK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_negctrls <- FSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_negctrls / FSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final

## add percentage of reads belonging to ASVs shared with corresponding isolation negative controls to table
MDK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_isolation_negctrls <- MDK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_isolation_negctrls / MDK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
PSK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_isolation_negctrls <- PSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_isolation_negctrls / PSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
MXK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_isolation_negctrls <- MXK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_isolation_negctrls / MXK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
FSK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_isolation_negctrls <- FSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_isolation_negctrls / FSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final

## add percentage of reads belonging to ASVs shared with corresponding library preparation negative control to table
MDK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_library_preparation_negctrls <- MDK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_library_preparation_negctrls / MDK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
PSK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_library_preparation_negctrls <- PSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_library_preparation_negctrls / PSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
MXK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_library_preparation_negctrls <- MXK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_library_preparation_negctrls / MXK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
FSK_milk_ASV_isol_Seq_negctrl_reads$perc_reads_shared_with_library_preparation_negctrls <- FSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_shared_with_library_preparation_negctrls / FSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final

## ensure the tables are sorted by sample_ID
MDK_milk_ASV_isol_Seq_negctrl_reads <- MDK_milk_ASV_isol_Seq_negctrl_reads[order(MDK_milk_ASV_isol_Seq_negctrl_reads$sample_ID),]
PSK_milk_ASV_isol_Seq_negctrl_reads <- PSK_milk_ASV_isol_Seq_negctrl_reads[order(PSK_milk_ASV_isol_Seq_negctrl_reads$sample_ID),]
MXK_milk_ASV_isol_Seq_negctrl_reads <- MXK_milk_ASV_isol_Seq_negctrl_reads[order(MXK_milk_ASV_isol_Seq_negctrl_reads$sample_ID),]
FSK_milk_ASV_isol_Seq_negctrl_reads <- FSK_milk_ASV_isol_Seq_negctrl_reads[order(FSK_milk_ASV_isol_Seq_negctrl_reads$sample_ID),]

## add column with pretty sample names to tables
MDK_milk_ASV_isol_Seq_negctrl_reads$sample <- c("MilkDNAKit_Milk-3_1", "MilkDNAKit_Milk-3_2", "MilkDNAKit_Milk-2_1", "MilkDNAKit_Milk-2_2", "MilkDNAKit_Milk-1_1", "MilkDNAKit_Milk-1_2")
PSK_milk_ASV_isol_Seq_negctrl_reads$sample <- c("PowerSoilProKit_Milk-3_1", "PowerSoilProKit_Milk-3_2", "PowerSoilProKit_Milk-3_3", "PowerSoilProKit_Milk-2_1", "PowerSoilProKit_Milk-2_2", "PowerSoilProKit_Milk-2_3", "PowerSoilProKit_Milk-1_1", "PowerSoilProKit_Milk-1_2", "PowerSoilProKit_Milk-1_3")
MXK_milk_ASV_isol_Seq_negctrl_reads$sample <- c("MagMAXKit_Milk-3_1", "MagMAXKit_Milk-3_2", "MagMAXKit_Milk-2_1", "MagMAXKit_Milk-2_2", "MagMAXKit_Milk-1_1", "MagMAXKit_Milk-1_2")
FSK_milk_ASV_isol_Seq_negctrl_reads$sample <- c("FastDNAStoolMiniKit_Milk-3_1", "FastDNAStoolMiniKit_Milk-3_2", "FastDNAStoolMiniKit_Milk-2_1", "FastDNAStoolMiniKit_Milk-2_2", "FastDNAStoolMiniKit_Milk-1_1", "FastDNAStoolMiniKit_Milk-1_2")

## resort rows and columns in tables
MDK_milk_ASV_isol_Seq_negctrl_reads <- MDK_milk_ASV_isol_Seq_negctrl_reads[order(MDK_milk_ASV_isol_Seq_negctrl_reads$sample), c(c(19,24,20,16:18,21:23), order(colnames(MDK_milk_ASV_isol_Seq_negctrl_reads)[1:15]))]
PSK_milk_ASV_isol_Seq_negctrl_reads <- PSK_milk_ASV_isol_Seq_negctrl_reads[order(PSK_milk_ASV_isol_Seq_negctrl_reads$sample), c(c(8,13,9,5:7,10:12), order(colnames(PSK_milk_ASV_isol_Seq_negctrl_reads)[1:4]))]
MXK_milk_ASV_isol_Seq_negctrl_reads <- MXK_milk_ASV_isol_Seq_negctrl_reads[order(MXK_milk_ASV_isol_Seq_negctrl_reads$sample), c(c(7,12,8,4:6,9:11), order(colnames(MXK_milk_ASV_isol_Seq_negctrl_reads)[1:3]))]
FSK_milk_ASV_isol_Seq_negctrl_reads <- FSK_milk_ASV_isol_Seq_negctrl_reads[order(FSK_milk_ASV_isol_Seq_negctrl_reads$sample), c(c(22,27,23,19:21,24:26), order(colnames(FSK_milk_ASV_isol_Seq_negctrl_reads)[1:18]))]

## calculate relative abundances
MDK_milk_ASV_isol_Seq_negctrl_reads[,10:ncol(MDK_milk_ASV_isol_Seq_negctrl_reads)] <- MDK_milk_ASV_isol_Seq_negctrl_reads[,10:ncol(MDK_milk_ASV_isol_Seq_negctrl_reads)] / MDK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
PSK_milk_ASV_isol_Seq_negctrl_reads[,10:ncol(PSK_milk_ASV_isol_Seq_negctrl_reads)] <- PSK_milk_ASV_isol_Seq_negctrl_reads[,10:ncol(PSK_milk_ASV_isol_Seq_negctrl_reads)] / PSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
MXK_milk_ASV_isol_Seq_negctrl_reads[,10:ncol(MXK_milk_ASV_isol_Seq_negctrl_reads)] <- MXK_milk_ASV_isol_Seq_negctrl_reads[,10:ncol(MXK_milk_ASV_isol_Seq_negctrl_reads)] / MXK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final
FSK_milk_ASV_isol_Seq_negctrl_reads[,10:ncol(FSK_milk_ASV_isol_Seq_negctrl_reads)] <- FSK_milk_ASV_isol_Seq_negctrl_reads[,10:ncol(FSK_milk_ASV_isol_Seq_negctrl_reads)] / FSK_milk_ASV_isol_Seq_negctrl_reads$n_reads_dada2_final

## change colnames for final read count
colnames(MDK_milk_ASV_isol_Seq_negctrl_reads)[3] <- "n_reads_total"
colnames(PSK_milk_ASV_isol_Seq_negctrl_reads)[3] <- "n_reads_total"
colnames(MXK_milk_ASV_isol_Seq_negctrl_reads)[3] <- "n_reads_total"
colnames(FSK_milk_ASV_isol_Seq_negctrl_reads)[3] <- "n_reads_total"

## save supplementary tables showing which bacteria / ASVs are detected in milk and their corresponding negative controls (1 table per DNA isolation method)
# write.table(MDK_milk_ASV_isol_Seq_negctrl_reads[,-1], file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Table_S4B.txt", sep="\t", row.names=F, quote=F)
# write.table(PSK_milk_ASV_isol_Seq_negctrl_reads[,-1], file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Table_S4C.txt", sep="\t", row.names=F, quote=F)
# write.table(MXK_milk_ASV_isol_Seq_negctrl_reads[,-1], file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Table_S4D.txt", sep="\t", row.names=F, quote=F)
# write.table(FSK_milk_ASV_isol_Seq_negctrl_reads[,-1], file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Table_S4A.txt", sep="\t", row.names=F, quote=F)

## combine info on number and percentage of reads belonging to ASVs detected in negative controls and add columns to original data frame
library(plyr)
milk_ASV_isol_Seq_negctrl_reads <- rbind.fill(MDK_milk_ASV_isol_Seq_negctrl_reads[,c(1,4:9)],
                                              PSK_milk_ASV_isol_Seq_negctrl_reads[,c(1,4:9)],
                                              MXK_milk_ASV_isol_Seq_negctrl_reads[,c(1,4:9)],
                                              FSK_milk_ASV_isol_Seq_negctrl_reads[,c(1,4:9)])
milk_meta <- left_join(milk_meta, milk_ASV_isol_Seq_negctrl_reads, by="sample_ID")


### MILK SAMPLES -> check which milk ASVs are uniquely detected in milk and not in any negative control

## identify ASVs unique to milk samples (i.e. not shared with the corresponding isolation and the library prep negative control)
MDK_milk_ASV_unique <- MDK_milk_ASV[,!(colnames(MDK_milk_ASV) %in% colnames(MDK_negctrl_ASV)) & !(colnames(MDK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] #141 unique ASVs
PSK_milk_ASV_unique <- PSK_milk_ASV[,!(colnames(PSK_milk_ASV) %in% colnames(PSK_negctrl_ASV)) & !(colnames(PSK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] #134 unique ASVs
MXK_milk_ASV_unique <- MXK_milk_ASV[,!(colnames(MXK_milk_ASV) %in% colnames(MXK_negctrl_ASV)) & !(colnames(MXK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] #176 unique ASVs
FSK_milk_ASV_unique <- FSK_milk_ASV[,!(colnames(FSK_milk_ASV) %in% colnames(FSK_negctrl_ASV)) & !(colnames(FSK_milk_ASV) %in% colnames(Seq_negctrl_ASV))] #173 unique ASVs

## add column showing sample ID
MDK_milk_ASV_unique$sample_ID <- rownames(MDK_milk_ASV_unique)
PSK_milk_ASV_unique$sample_ID <- rownames(PSK_milk_ASV_unique)
MXK_milk_ASV_unique$sample_ID <- rownames(MXK_milk_ASV_unique)
FSK_milk_ASV_unique$sample_ID <- rownames(FSK_milk_ASV_unique)

MDK_milk_ASV_unique <- MDK_milk_ASV_unique[,c(142,1:141)] #resort columns in first data frame so the sample_ID column shows first

## combine in 1 data frame
milk_ASV_unique <- rbind.fill(MDK_milk_ASV_unique, PSK_milk_ASV_unique, MXK_milk_ASV_unique, FSK_milk_ASV_unique)
milk_ASV_unique <- milk_ASV_unique[order(milk_ASV_unique$sample_ID),] #sort by sample_ID
rownames(milk_ASV_unique) <- milk_ASV_unique$sample_ID
milk_ASV_unique <- milk_ASV_unique[,-1] #remove sample_ID column
milk_ASV_unique[is.na(milk_ASV_unique)] <- 0 #if entries were filled up that was done with NA, change them to 0 as the entries should show 0 reads for that column

## subset unique-to-milk data by taxonomic level
# note that these data frames still contain ASVs, which are unassigned on some taxonomic levels, incl. on kingdom and phylum level
milk_genus_unique <- subset.by.taxlevel(inputdata = milk_ASV_unique, pattern = "_s__.*") #27x159 genera (8 genera less)
#milk_ASV_unique #27x393 ASVs


### ===== 5. MILK SAMPLES ONLY: CALCULATE RICHNESS AFTER EXCLUDING CONTAMINATION ===== ###

#for each milk samples, calculate richness for genus level and for ASV level - only taking into account bacteria that could be classified on the given level

### MILK

## calculate richness for different taxonomic levels
milk_genus_unique_richness <- as.data.frame(specnumber(milk_genus_unique[,grep("g__NA", colnames(milk_genus_unique), invert=T)]))
milk_ASV_unique_richness <- as.data.frame(specnumber(milk_ASV_unique[,grep("ASV__NA", colnames(milk_ASV_unique), invert=T)]))

## combine in 1 data frame
milk_unique_richness <- as.data.frame(cbind(milk_genus_unique_richness, milk_ASV_unique_richness))
milk_unique_richness$sample_ID <- rownames(milk_unique_richness)
milk_unique_richness <- milk_unique_richness[order(milk_unique_richness$sample_ID),]
milk_unique_richness$DNA_isolation_method <- c(rep("MilkDNAKit",6), rep("PowerSoilProKit",9), rep("MagMAXKit",6), rep("FastStoolMiniKit",6))
colnames(milk_unique_richness)[1:2] <- c("milk_unique_genera_richness", "milk_unique_ASV_richness")
milk_unique_richness <- milk_unique_richness[,c(3:4,1:2)]
milk_unique_richness

# add unique richness info to original data frames
rownames(milk_meta) <- milk_meta$sample_ID
milk_meta <- merge(milk_meta, milk_unique_richness[,-c(1:2)], by="row.names")
milk_meta <- milk_meta[,-1]


### ===== 6. MOCK AND MILK SAMPLES ONLY: CALCULATE PERCENTAGE OF READS CLASSIFIED ON SPECIES LEVEL ===== ###

### MOCK SAMPLES
## save data that could be classified on the species level separately
mock_species_classified <- as.data.frame(mock_species[,grep("s__NA", colnames(mock_species), invert=T)])

## calculate number of reads that could be classified on the species level
mock_species_classified$n_reads_classified_on_species_level <- rowSums(mock_species_classified)

## add info to mock_meta data frame
mock_species_classified$sample_ID <- rownames(mock_species_classified)
mock_meta <- left_join(mock_meta, mock_species_classified[,c(4:5)], by="sample_ID")

## calculate percentage of reads that could be classified on the species level
mock_meta$perc_reads_classified_on_species_level <- mock_meta$n_reads_classified_on_species_level / mock_meta$n_reads_dada2_final


### ===== 7. TABLE 1A: SUMMARY STATISTICS FOR MOCK SAMPLES ===== ###

# split mock_meta data frame by DNA isolation method
mock_meta_MDK <- mock_meta[mock_meta$DNA_isolation_method=="MilkDNAKit" & !is.na(mock_meta$DNA_isolation_method),]
mock_meta_PSK <- mock_meta[mock_meta$DNA_isolation_method=="PowerSoilProKit" & !is.na(mock_meta$DNA_isolation_method),]
mock_meta_MXK <- mock_meta[mock_meta$DNA_isolation_method=="MagMAXKit" & !is.na(mock_meta$DNA_isolation_method),]
mock_meta_FSK <- mock_meta[mock_meta$DNA_isolation_method=="FastStoolMiniKit" & !is.na(mock_meta$DNA_isolation_method),]
mock_meta_DNA_Mock <- mock_meta[is.na(mock_meta$DNA_isolation_method),]

# function to retrieve values
getSumStats_mock <- function(inputdata){
  SumStats_mock <- c(nrow(inputdata),
                      inputdata$DNA_yield_ng,
                      inputdata$n_reads_dada2_final,
                      inputdata$ASV_richness,
                      inputdata$genera_richness,
                      inputdata$perc_reads_classified_on_species_level*100)
  return(SumStats_mock)
}

SumStats_mock <- as.data.frame(cbind(getSumStats_mock(mock_meta_MDK),
                                     getSumStats_mock(mock_meta_PSK),
                                     getSumStats_mock(mock_meta_MXK),
                                     getSumStats_mock(mock_meta_FSK),
                                     getSumStats_mock(mock_meta_DNA_Mock)))
SumStats_mock[7,] <- c(rep("Bacterial mock community",4), "DNA mock community")
SumStats_mock[8,] <- c("Milk Bacterial DNA Isolation Kit","PowerSoil Pro Kit","MagMAX Kit","Fast Stool DNA Mini Kit","NA")
SumStats_mock$info <- c("n (samples)", "DNA yield (ng)", "n (reads)", "n (ASVs)", "n (classified genera)", "% (reads classified on species level)", "Sample", "DNA isolation method")
SumStats_mock <- SumStats_mock[c(7:8,1:6), c(6,4,1,2,3,5)]
SumStats_mock
# write.table(SumStats_mock, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Table_1A.txt", sep="\t", row.names=F, quote=F)


### ===== 8. TABLE 1B: SUMMARY STATISTICS FOR NEGATIVE CONTROLS ===== ###

# split negctrl_meta data frame by DNA isolation method
negctrl_meta_MDK <- negctrl_meta[negctrl_meta$DNA_isolation_method=="MilkDNAKit" & !is.na(negctrl_meta$DNA_isolation_method),]
negctrl_meta_PSK <- negctrl_meta[negctrl_meta$DNA_isolation_method=="PowerSoilProKit" & !is.na(negctrl_meta$DNA_isolation_method),]
negctrl_meta_MXK <- negctrl_meta[negctrl_meta$DNA_isolation_method=="MagMAXKit" & !is.na(negctrl_meta$DNA_isolation_method),]
negctrl_meta_FSK <- negctrl_meta[negctrl_meta$DNA_isolation_method=="FastStoolMiniKit" & !is.na(negctrl_meta$DNA_isolation_method),]
negctrl_meta_Seq <- negctrl_meta[is.na(negctrl_meta$DNA_isolation_method),]

# function to retrieve median (range)
getSumStats_negctrl <- function(inputdata){
  SumStats_negctrl <- c(nrow(inputdata),
                      paste0(median(inputdata$DNA_yield_ng), " (", min(inputdata$DNA_yield_ng), " - ", max(inputdata$DNA_yield_ng), ")"),
                      paste0(median(inputdata$n_reads_dada2_final), " (", min(inputdata$n_reads_dada2_final), " - ", max(inputdata$n_reads_dada2_final), ")"),
                      paste0(median(inputdata$ASV_richness, na.rm=T), " (", min(inputdata$ASV_richness, na.rm=T), " - ", max(inputdata$ASV_richness, na.rm=T), ")"), #note: remove the negctrl with 0 reads here
                      paste0(median(inputdata$genera_richness, na.rm=T), " (", min(inputdata$genera_richness, na.rm=T), " - ", max(inputdata$genera_richness, na.rm=T), ")")) #note: remove the negctrl with 0 reads here
  return(SumStats_negctrl)
}

SumStats_negctrl <- as.data.frame(cbind(getSumStats_negctrl(negctrl_meta_MDK),
                                        getSumStats_negctrl(negctrl_meta_PSK),
                                        getSumStats_negctrl(negctrl_meta_MXK),
                                        getSumStats_negctrl(negctrl_meta_FSK),
                                        getSumStats_negctrl(negctrl_meta_Seq)))
for (i in 1:ncol(SumStats_negctrl)){SumStats_negctrl[,i] <- as.character(as.factor(SumStats_negctrl[,i]))}
SumStats_negctrl[6,] <- c(rep("DNA isolation negative controls",4), "Library preparation negative control")
SumStats_negctrl[7,] <- c("Milk Bacterial DNA Isolation Kit","PowerSoil Pro Kit","MagMAX Kit","Fast Stool DNA Mini Kit","NA")
SumStats_negctrl$info <- c("n (samples)", "DNA yield (ng)", "n (reads)", "n (ASVs)", "n (classified genera)", "Samples", "DNA isolation method")
SumStats_negctrl <- SumStats_negctrl[c(6,7,1:5), c(6,4,1,2,3,5)]
SumStats_negctrl
# write.table(SumStats_negctrl, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Table_1B.txt", sep="\t", row.names=F, quote=F)


### ===== 9. TABLE 1C: SUMMARY STATISTICS FOR MILK SAMPLES ===== ###

## split milk_meta data frame by DNA isolation method
milk_meta_MDK <- milk_meta[milk_meta$DNA_isolation_method=="MilkDNAKit",]
milk_meta_PSK <- milk_meta[milk_meta$DNA_isolation_method=="PowerSoilProKit",]
milk_meta_MXK <- milk_meta[milk_meta$DNA_isolation_method=="MagMAXKit",]
milk_meta_FSK <- milk_meta[milk_meta$DNA_isolation_method=="FastStoolMiniKit",]

## function to retrieve median (range)
getSumStats_milk <- function(inputdata){
  SumStats_milk <- c(nrow(inputdata),
                     paste0(median(inputdata$DNA_yield_ng), " (", min(inputdata$DNA_yield_ng), " - ", max(inputdata$DNA_yield_ng), ")"),
                     paste0(median(inputdata$n_reads_dada2_final), " (", min(inputdata$n_reads_dada2_final), " - ", max(inputdata$n_reads_dada2_final), ")"),
                     paste0(median(inputdata$ASV_richness), " (", min(inputdata$ASV_richness), " - ", max(inputdata$ASV_richness), ")"),
                     paste0(median(inputdata$genera_richness), " (", min(inputdata$genera_richness), " - ", max(inputdata$genera_richness), ")"),
                     paste0(median(inputdata$milk_unique_genera_richness), " (", min(inputdata$milk_unique_genera_richness), " - ", max(inputdata$milk_unique_genera_richness), ")"),
                     paste0(median(inputdata$perc_reads_shared_with_negctrls)*100, "% (", min(inputdata$perc_reads_shared_with_negctrls)*100, "% - ", max(inputdata$perc_reads_shared_with_negctrls)*100, "%)"))
  
  names(SumStats_milk) <- c("n (samples)", "DNA yield (ng)", "n (reads)", "n (ASVs)", "n (classified genera)", "n (classified genera unique to milk)", "% (reads shared with negative controls)")
  
  return(SumStats_milk)
}

SumStats_milk <- as.data.frame(cbind(getSumStats_milk(milk_meta_MDK),
                                        getSumStats_milk(milk_meta_PSK),
                                        getSumStats_milk(milk_meta_MXK),
                                        getSumStats_milk(milk_meta_FSK)))
SumStats_milk$info <- rownames(SumStats_milk)

# function to get p and FDR
getStats_milk <- function(inputdata, columns_of_interest){
  p <- c()
  
  for (i in columns_of_interest){
    print(colnames(inputdata)[i])
    kw <- kruskal.test(inputdata[,i] ~ DNA_isolation_method, data=inputdata)
    p <- c(p, kw$p.value)
  }
  
  FDR <- p.adjust(p, method = "BH")
  results <- data.frame(p=round(p,4), FDR=round(FDR,3))
  rownames(results)<- c("DNA yield (ng)", "n (reads)", "n (ASVs)", "n (classified genera)", "n (classified genera unique to milk)", "% (reads shared with negative controls)")
  
  return(results)
}

Stats_milk <- getStats_milk(milk_meta, columns_of_interest = c(16,18,20,19,27,24))
Stats_milk$info <- rownames(Stats_milk)

## add statistics to milk summary stats table
SumStats_milk2 <- left_join(SumStats_milk[,c(5,1:4)], Stats_milk, by="info")

## clean up milk summary stats table
for (i in 1:ncol(SumStats_milk2)){SumStats_milk2[,i] <- as.character(as.factor(SumStats_milk2[,i]))}

SumStats_milk2[8,] <- c("Samples", "3 milk samples in duplicate", "3 milk samples in triplicate", "3 milk samples in duplicate", "3 milk samples in duplicate", NA, NA)
SumStats_milk2[9,] <- c("DNA isolation method", "Milk Bacterial DNA Isolation Kit","PowerSoil Pro Kit","MagMAX Kit","Fast Stool DNA Mini Kit", NA, NA)

SumStats_milk2 <- SumStats_milk2[c(8,9,1:7),c(1,5,2:4,6:7)]
SumStats_milk2

## add 16S rRNA gene copies per 1ml milk (only for PSK-isolated milk samples)
summary(milk_meta[milk_meta$DNA_isolation_method=="PowerSoilProKit" & !is.na(milk_meta$qPCR_16S_rRNA_gene_copies_per_1ml_milk_input_log10),"qPCR_16S_rRNA_gene_copies_per_1ml_milk_input_log10"])
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 4.813   4.891   4.952   5.017   5.213   5.218
SumStats_milk2[10,] <- c("n (16S rRNA gene copies per 1 ml milk)", NA, NA, "10^5.0 (10^4.8 - 10^5.2)", NA, NA, NA)
SumStats_milk2

## save final Table 1C
# write.table(SumStats_milk2, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Table_1C.txt", sep="\t", row.names=F, quote=F)


## run Dunn's posthoc tests for parameters with FDR<0.05
library(FSA)

# Dunn's posthoc test for DNA yield
DT1 <- dunnTest(DNA_yield_ng ~ DNA_isolation_method, data=milk_meta, method="bh")
dunnTestRes_DNA_yield <- as.data.frame(DT1$res)
# write.table(dunnTestRes_DNA_yield, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Dunns_test_results_DNA_yield.txt", sep="\t", row.names=F, quote=F)

# Dunn's posthoc test for read count
DT2 <- dunnTest(n_reads_dada2_final ~ DNA_isolation_method, data=milk_meta, method="bh")
dunnTestRes_reads <- as.data.frame(DT2$res)
# write.table(dunnTestRes_reads, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Dunns_test_results_reads.txt", sep="\t", row.names=F, quote=F)

# Dunn's posthoc test for ASV richness
DT3 <- dunnTest(ASV_richness ~ DNA_isolation_method, data=milk_meta, method="bh")
dunnTestRes_ASV_richness <- as.data.frame(DT3$res)
# write.table(dunnTestRes_ASV_richness, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Dunns_test_results_ASV_richness.txt", sep="\t", row.names=F, quote=F)


### ===== 10. COMPARISON OF READ COUNT BETWEEN MILK AND NEGATIVE CONTROLS ===== ###

## create dataset for test
reads <- data.frame(rbind(negctrl_meta[,c(13,14,18)], milk_meta[,c(13,14,18)]))

## get summary statistics
reads_negctrl <- summary(reads[reads$sample_type=="Negative_control",3])
reads_milk <- summary(reads[reads$sample_type=="Milk",3])

## wilcoxon test
WT <- wilcox.test(reads[reads$sample_type=="Negative_control",3],
                  reads[reads$sample_type=="Milk",3],
                  alternative="two.sided", paired=F, exact=T)

## create and save output table
reads_table <- data.frame(median_n_reads_negctrl=round(median(reads[reads$sample_type=="Negative_control",3]),0),
                          median_n_reads_milk=round(median(reads[reads$sample_type=="Milk",3]),0),
                          p=WT$p.value)
# write.table(reads_table, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Wilcoxon_test_results_read_count_milk_vs_negctrls.txt", sep="\t", row.names=F, quote=F)


### ===== 11. CORRELATION OF READ COUNT AND DNA YIELD WITH PERCENTAGE OF CONTAMINATION IN MILK SAMPLES ===== ###

## check data distribution
hist(milk_meta$n_reads_dada2_final, breaks=20)
hist(milk_meta$perc_reads_shared_with_negctrls, breaks=20)
hist(milk_meta$perc_reads_shared_with_isolation_negctrls, breaks=20)
hist(milk_meta$perc_reads_shared_with_library_preparation_negctrls, breaks=20)
hist(milk_meta$DNA_yield_ng, breaks=20)
# data is left-skewed -> use non-parametric tests

## correlation plot for read count~perc of contamination
# all samples from all DNA isolation methods together
library(ggpubr)
ggscatter(milk_meta, x="n_reads_dada2_final", y="perc_reads_shared_with_negctrls", add="reg.line", cor.method="spearman", title="Spearman correlation")+stat_cor(method="spearman", label.x=30) #rho=-0.0061, p=0.98

## correlation plot for DNA yield~perc of contamination
# all samples from all DNA isolation methods together
ggscatter(milk_meta, x="DNA_yield_ng", y="perc_reads_shared_with_negctrls", add="reg.line", cor.method="spearman", title="Spearman correlation")+stat_cor(method="spearman", label.x=30) #rho=-0.56, p=0.0026

## save resutls in data frame and run FDR correction
corresults <- data.frame(x=c("n_reads_dada2_final", "DNA_yield_ng"),
                         y=c("perc_reads_shared_with_negctrls","perc_reads_shared_with_negctrls"),
                         rho=c(-0.0061,-0.56),
                         p=c(0.98,0.0026))
corresults$FDR <- p.adjust(corresults$p, method="BH")
corresults
# write.table(reads_table, file="MY_PATH_TO_DNA_ISOLATION_COMPARISON_RESULTS/Correlation_results.txt", sep="\t", row.names=F, quote=F)














