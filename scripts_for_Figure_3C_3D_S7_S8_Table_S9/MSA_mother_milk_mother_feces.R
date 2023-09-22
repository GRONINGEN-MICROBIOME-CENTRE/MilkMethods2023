################################################################################
##### Per-genus multiple sequence alignments of ASVs in maternal milk and maternal fecal samples
### Author(s): Asier Fernández / Johanne E. Spreckels
### Last updated: 20th September, 2023
################################################################################

#****************
# Load libraries
#****************
library(stringr)
library(seqinr)
library(msa)

#********************
# Defining functions
#********************
#A. MSA_result: generate a MSA of the ASVs of selected genera using ClustalW 
# It outputs a FASTA file with ASVs sequences, a FASTA file with the alignment and a PDF with the alignment (including logos)
# output_dir : full path to directory where one folder will be stored for each genus with its results
# description: information about sequencing technology and pairs of samples to compare
MSA_result <- function(ASVs, names_ASVs, output_dir, description) {
  for (i in 1:length(names_ASVs)) {
    genus <- names_ASVs [i] #name of the genus
    dir <- paste0(output_dir ,genus)
    dir.create(dir) #Create a directory for the genus
    setwd (dir)  #Move to the directory
    genus_ASVs <- ASVs[grep (genus, names(ASVs))] #select ASVs of the specific genus
    filename <- paste0(paste(genus,description, sep = "_"), ".fa")
    write.fasta(genus_ASVs, names(genus_ASVs), filename) #generate a FASTA formatted file
    MSA <-  msa(filename, method = "ClustalW", type = "DNA", order = "input") #do the MSA
    ASV_names_for_plot <- gsub("\\..*", "", names(genus_ASVs)) #get names of the ASVs sequences to be displayed later in the alignment and the tree
    names(MSA@unmasked) <- ASV_names_for_plot
    try({msaPrettyPrint(MSA, output="pdf", file = paste(genus, description, "MSA.pdf", sep= "_"), #get a PDF with the alignment
                   showNames ="left", shadingMode="similar",
                   shadingColors="blues", showLogo="top", logoColors="rasmol",
                   showLegend=T, askForOverwrite=FALSE)})
    msa_align <- msaConvert(MSA, type="seqinr::alignment") #convert the alignment for later processing with seqinr
    msa_align$nam <- ASV_names_for_plot #add names for the plot
    filename_aln <- paste0(paste(genus,description, sep = "_"), "_aln.fa")
    write.fasta(as.list(msa_align$seq), names(genus_ASVs), filename_aln) #write alignment fasta (as ALN file is empty)
    }
}

# Set working directory
setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/")

#Read metadata
metadata <- read.delim("Data_sequencing_method_comparison_n136.txt") 


######################################################################
# Genus-level MSAs of ASVs in maternal milk-fecal samples
######################################################################

#################################################
#A. 16S V3-V4 - Mother Milk vs Mother Feces 
#################################################
setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/16S/Mother_Milk-Mother_Feces/")

#*****************
# A1.Get metadata
#*****************
metadata_16S <- metadata[grepl("pilot", metadata$description) &
                           metadata$sample_origin == "Mother" &
                           metadata$sequencing_method == "16S_V3V4" ,]

# Remove the pair with 0 reads (in a maternal fecal sample)
pair_to_exclude <- metadata_16S$pair_ID[which(metadata_16S$n_reads_total == 0)]
metadata_16S <- metadata_16S[!(metadata_16S$pair_ID == pair_to_exclude),]


#***********************************************************************
# A2.Get the genera (and their ASVs) present in at least 2 paired samples
#***********************************************************************
# For each sample select those ASVs that are present
# Use the sequencing technology +  sample_ID + lowest-level taxonomic assigment as ASV name
# Concatenate all sequences in a list
sample_ASVs_16S <- metadata_16S[,22:ncol(metadata_16S)]
rownames(sample_ASVs_16S) <- paste(metadata_16S$sequencing_method,metadata_16S$sample_origin, metadata_16S$sample_type, metadata_16S$pair_ID, sep = "_")
ASV_final_list_16S <- list()

for (i in 1:nrow(sample_ASVs_16S)) {
  sample_name <- rownames(sample_ASVs_16S)[i]
  ASV <- gsub(".*_", "",  colnames(sample_ASVs_16S)[which(sample_ASVs_16S[i,] != 0)]) #sequence
  name_full <- gsub("_ASV.*", "\\1", colnames(sample_ASVs_16S)[which(sample_ASVs_16S[i,] != 0)]) #full-name
  name <- gsub("_s__.*", "", name_full) #remove everithing after genus level
  name <- gsub(".*g__", "g__", name) #remove everything before genus level (get genus name)
  indexes_to_update <- which(str_sub(name, - 2, - 2) == "_") #get indexes of genera names that need to be updated (extra undescore in the end)
  name[indexes_to_update] <- sub("_[^_]+$", "", name[indexes_to_update])  #update names    
  if (any(grepl("__NA",name_full))) { #if genus is missing get the lowest-level taxonomy
    NAs <- which(name == "g__NA")
    name [NAs] <-  gsub("..__NA.*", "", name_full[which(name == "g__NA")])
    name [NAs] <-  gsub(".*_", "", name [NAs])
  }
  ASV_list <- ASV
  names(ASV_list) <- paste(rep(sample_name, length(name)), name, sep = ".")
  ASV_final_list_16S <- append(ASV_final_list_16S, ASV_list) #final list of ASVs with names and sequences
}

# Check ASV names in case there are genera that need to be renamed: names(ASV_final_list_16S)

# Create a dataframe with the counts of each genus for each participant 
names_genera_16S <- gsub(".*g__", "", grep("g__", names(ASV_final_list_16S), value = T))
names_genera_16S <- unique(names_genera_16S) # 126 unique genera
genera_counts_16S <- data.frame(matrix(ncol = length(rownames(sample_ASVs_16S)) + 1, nrow = length(names_genera_16S)))
colnames(genera_counts_16S) <- c("genera", rownames(sample_ASVs_16S))
genera_counts_16S$genera <- names_genera_16S

for (i in 2:ncol(genera_counts_16S)) {
  microbes <- names(ASV_final_list_16S)[grep(paste0(colnames(genera_counts_16S)[i],"\\."), names(ASV_final_list_16S))]
  genera <- gsub(".*g__", "", grep("g__", microbes, value = T))
  counts <- data.frame(table(genera))
  genera_counts_16S[i] <- counts$Freq[match(genera_counts_16S$genera,counts$genera)]
}
genera_counts_16S[is.na(genera_counts_16S)] <- 0 #set NAs to 0

# Create a list with pair_IDs as elements and that each element has the names of the genera present in both sample types
genera_in_pairs_16S <- list()
IDs <- unique(metadata_16S$pair_ID)
for (i in IDs) {
  index <- which(IDs == i)
  pattern <- paste0("_",i, "$")
  pair_counts <- genera_counts_16S[, grep(pattern, colnames(genera_counts_16S))]
  genera_in_pairs_16S [[index]] <- genera_counts_16S$genera [as.numeric(rownames(pair_counts[apply(pair_counts!=0, 1, all),]))] #select genera detected in at least 1 pair 
}

genera_in_pairs_16S <- genera_in_pairs_16S[!sapply(genera_in_pairs_16S,is.null)] # removes null elements
names(genera_in_pairs_16S) <- IDs

# Get names of the shared genera present in 2 or more pairs
shared_genera_16S <- names(which(table(unlist(genera_in_pairs_16S)) > 1))

# Get all ASVs of selected genera
names(ASV_final_list_16S) <- make.unique(names(ASV_final_list_16S), sep = '_') # distinguish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
matching_ASVs_16S <- ASV_final_list_16S[unique (grep(paste(shared_genera_16S,collapse="|"), # get all ASVs belonging to shared genera
                                             names(ASV_final_list_16S)))]

#***********************************************************************
# A3.Get the MSA for each genus that is detected among paired individuals
#***********************************************************************
MSA_result(matching_ASVs_16S,shared_genera_16S, "~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/16S/Mother_Milk-Mother_Feces/", "16S_mother_milk_mother_feces")

#################################################
#B. 16S_23S - Mother Milk vs Mother Feces
#################################################

setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/16S23S/Mother_Milk-Mother_Feces//")

#*****************
# B1.Get metadata
#*****************
metadata_16S23S <- metadata[grepl("pilot", metadata$description) &
                           metadata$sample_origin == "Mother" &
                           metadata$sequencing_method == "16S-ITS-23S" ,]

#***********************************************************************
# B2.Get the genera (and their ASVs) present in at least 2 paired samples
#***********************************************************************
# For each sample select those ASVs that are present
# Use the sequencing technology +  sample_ID + lowest-level taxonomic assigment as ASV name
# Concatenate all sequences in a list
sample_ASVs_16S23S <- metadata_16S23S[,22:ncol(metadata_16S23S)]
rownames(sample_ASVs_16S23S) <- paste(metadata_16S23S$sequencing_method,metadata_16S23S$sample_origin, metadata_16S23S$sample_type, metadata_16S23S$pair_ID, sep = "_")
ASV_final_list_16S23S <- list()

for (i in 1:nrow(sample_ASVs_16S23S)) {
  sample_name <- rownames(sample_ASVs_16S23S)[i]
  ASV <- gsub(".*_", "",  colnames(sample_ASVs_16S23S)[which(sample_ASVs_16S23S[i,] != 0)]) #sequence
  name_full <- gsub("_ASV.*", "\\1", colnames(sample_ASVs_16S23S)[which(sample_ASVs_16S23S[i,] != 0)]) #full-name
  name <- gsub("_s__.*", "", name_full) #remove everything after genus level
  name <- gsub(".*g__", "g__", name) #remove everything before genus level (get genus name)
  indexes_to_update <- which(str_sub(name, - 2, - 2) == "_") #get indexes of genera names that need to be updated (extra undescore in the end)
  name[indexes_to_update] <- sub("_[^_]+$", "", name[indexes_to_update])  #update names    
  if (any(grepl("__NA",name_full))) { #if genus is missing get the lowest-level taxonomy
    NAs <- which(name == "g__NA")
    name [NAs] <-  gsub("..__NA.*", "", name_full[which(name == "g__NA")])
    name [NAs] <-  gsub(".*_", "", name [NAs])
  }
  ASV_list <- ASV
  names(ASV_list) <- paste(rep(sample_name, length(name)), name, sep = ".")
  ASV_final_list_16S23S <- append(ASV_final_list_16S23S, ASV_list) #final list of ASVs with names and sequences
}

#Check ASV names in case there are genera that need to be renamed: names(ASV_final_list_16S23S)

#Create a dataframe with the counts of each genus for each participant 
names_genera_16S23S <- gsub(".*g__", "", grep("g__", names(ASV_final_list_16S23S), value = T))
names_genera_16S23S <- unique(names_genera_16S23S) #57
genera_counts_16S23S <- data.frame(matrix(ncol = length(rownames(sample_ASVs_16S23S)) + 1, nrow = length(names_genera_16S23S)))
colnames(genera_counts_16S23S) <- c("genera", rownames(sample_ASVs_16S23S))
genera_counts_16S23S$genera <- names_genera_16S23S

for (i in 2:ncol(genera_counts_16S23S)) {
  microbes <- names(ASV_final_list_16S23S)[grep(paste0(colnames(genera_counts_16S23S)[i],"\\."), names(ASV_final_list_16S23S))]
  genera <- gsub(".*g__", "", grep("g__", microbes, value = T))
  counts <- data.frame(table(genera))
  genera_counts_16S23S[i] <- counts$Freq[match(genera_counts_16S23S$genera,counts$genera)]
}

genera_counts_16S23S[is.na(genera_counts_16S23S)] <- 0 #set NAs to 0
cbind(genera_counts_16S23S$genera, rowSums(genera_counts_16S23S[2:ncol(genera_counts_16S23S)])) #check how many times each genera is detected

# Create a list with pair_IDs as elements and that each element has the names of the genera present in both sample types
genera_in_pairs_16S23S <- list()
IDs <- unique(metadata_16S23S$pair_ID)
for (i in IDs) {
  index <- which(IDs == i)
  pattern <- paste0("_",i, "$")
  pair_counts <- genera_counts_16S23S[, grep(pattern, colnames(genera_counts_16S23S))]
  genera_in_pairs_16S23S [[index]] <- genera_counts_16S23S$genera [as.numeric(rownames(pair_counts[apply(pair_counts!=0, 1, all),]))] #select genera detected in at least 1 pair of oral-feces
}

genera_in_pairs_16S23S <- genera_in_pairs_16S23S[!sapply(genera_in_pairs_16S23S,is.null)] # removes null elements
names(genera_in_pairs_16S23S) <- IDs

# Get names of the shared genera present in 2 or more pairs
shared_genera_16S23S <- names(which(table(unlist(genera_in_pairs_16S23S)) > 1))

# Get all ASVs of selected genera
names(ASV_final_list_16S23S) <- make.unique(names(ASV_final_list_16S23S), sep = '_') #distingish between the names of different ASVs with same taxonomy (adding _1, _2, etc)
matching_ASVs_16S23S <- ASV_final_list_16S23S[unique(grep(paste(shared_genera_16S23S,collapse="|"), # get all ASVs belonging to shared genera
                                             names(ASV_final_list_16S23S)))]

#***********************************************************************
# B3.Get the MSA for each genus that is detected among paired individuals
#***********************************************************************
MSA_result(matching_ASVs_16S23S, shared_genera_16S23S, "~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/16S23S/Mother_Milk-Mother_Feces/", "16S23S_mother_milk_mother_feces")


########################
#C. Save results 
########################
setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/")

# 16S V3-V4 results
write.table(genera_counts_16S,"16S/Mother_Milk-Mother_Feces/Genera_counts.txt", sep = "\t", row.names = T, quote = FALSE) 
saveRDS(genera_in_pairs_16S, "16S/Mother_Milk-Mother_Feces/Genera_counts_pairs.rds")
saveRDS(matching_ASVs_16S, "16S/Mother_Milk-Mother_Feces/ASVs_shared_genera.rds")

# 16S-ITS-23S results
write.table(genera_counts_16S23S,"16S23S/Mother_Milk-Mother_Feces/Genera_counts.txt", sep = "\t", row.names = T, quote = FALSE) 
saveRDS(genera_in_pairs_16S23S, "16S23S/Mother_Milk-Mother_Feces/Genera_counts_pairs.rds")
saveRDS(matching_ASVs_16S23S, "16S23S/Mother_Milk-Mother_Feces/ASVs_shared_genera.rds")

