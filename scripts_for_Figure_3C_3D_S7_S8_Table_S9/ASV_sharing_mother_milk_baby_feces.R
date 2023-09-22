################################################################################
##### ASVs transmission analysis between maternal milk and infant fecal samples
### Author(s): Asier Fernández / Johanne E. Spreckels
### Last updated: 20th September, 2023
################################################################################

#********************
# Defining functions
#********************
# Function to estimate the number of related and unrelated distances == or != 0.
# It also estimates the number of pairs in which PSE are detected
# SNP_distances: Dataframe with SNP distances of the ASVs from the shared genera 
# ASVs: list of sequences of ASVs belonging to the specific genus (e.g. matching_ASVs_genus_16S23S)
exact_PSE_analysis <- function(SNP_distances, ASVs) {
  #Define vectors to store results
  counts_related_below <- 0
  counts_related_above <-0
  counts_unrelated_below<-0
  counts_unrelated_above <-0
  pairs_with_PSE<-c()
  unrelated_comp_with_PSE <- list()
  unrelated_comp_without_PSE <- list()
  matching_ASVs_genus <- ASVs[names(ASVs) %in% rownames(SNP_distances)] #select the ASVs from the specific genus
  # Iterate over the distance matrix
  for (i in 1:nrow(SNP_distances)) {
    for (j in 1:ncol(SNP_distances)) {
      if(!is.na(SNP_distances[i,j])) {
        ASV1_fullname <- rownames(SNP_distances)[i] #get full ASV same
        ASV1_sampleorigin <- gsub(".*_","", gsub("(.*)_\\w+", "\\1", gsub("(.*)_\\w+", "\\1", gsub("\\..*","",  ASV1_fullname)))) #will be mother/baby, except for oral samples that will be oral (we know its a baby then)
        ASV1_pair <- gsub(".*_","", gsub("\\..*","",  ASV1_fullname)) #i #get pair_ID
        ASV1_seq_length <- str_length(matching_ASVs_genus[ASV1_fullname]) #get sequence length
        ASV2_fullname <- colnames(SNP_distances) [j] # get full ASV name
        ASV2_pair <- gsub(".*_","",gsub("\\..*","",  ASV2_fullname)) # get pair ID
        ASV2_sampleorigin <- gsub(".*_","", gsub("(.*)_\\w+", "\\1", gsub("(.*)_\\w+", "\\1", gsub("\\..*","",  ASV2_fullname)))) 
        ASV2_seq_length <- str_length(matching_ASVs_genus[ASV2_fullname]) #get sequence length
        mean_seq_length <- mean(c(ASV1_seq_length, ASV2_seq_length)) #estimate the mean length of both ASVs
        if (SNP_distances[i,j] == 0) { #if the distance is below our threshold
          if (ASV1_sampleorigin != ASV2_sampleorigin) { #exclude comparisons from ASVs from the same sample origin
            if (ASV1_fullname != ASV2_fullname) { # exclude the same ASVs (obviously their distance is 0)
              if (ASV1_pair == ASV2_pair) { # if they are related (same pair_ID)
                pairs_with_PSE<-c(pairs_with_PSE, ASV1_pair) #add name of the pair
                counts_related_below <- counts_related_below + 1 #add 1 count to the related group
              } else { # if they are unrelated (different pair_ID)
                counts_unrelated_below <- counts_unrelated_below + 1 #add 1 count to the unrelated group
                unrelated_comp_with_PSE <-c(unrelated_comp_with_PSE, list(paste(paste(ASV1_sampleorigin,ASV1_pair, sep ="_"), 
                                                                                paste(ASV2_sampleorigin,ASV2_pair, sep ="_") , sep="_"))) #keep the unrelated pairs with PSE
                unrelated_comp_with_PSE <-c(unrelated_comp_with_PSE, list(paste(paste(ASV1_sampleorigin,ASV1_pair, sep ="_"), 
                                                                                paste(ASV2_sampleorigin,ASV2_pair, sep ="_") , sep="_"))) # all combinations
              }
            }
          }
        } else { #if the distance is above our threshold
          if (ASV1_sampleorigin != ASV2_sampleorigin) { #exclude comparisons from ASVs from the same sample origin
            if (ASV1_fullname != ASV2_fullname) { # exclude the same ASVs (obviously their distance is 0)
              if (ASV1_pair == ASV2_pair) { # if they are related (same pair_ID)
                counts_related_above <- counts_related_above + 1 #add 1 count to the related group
              } else { # if they are unrelated (different pair_ID)
                counts_unrelated_above <- counts_unrelated_above + 1 #add 1 count to the unrelated group
                unrelated_comp_without_PSE <-c(unrelated_comp_without_PSE, list(paste(paste(ASV1_sampleorigin,ASV1_pair, sep ="_"), 
                                                                                      paste(ASV2_sampleorigin,ASV2_pair, sep ="_") , sep="_"))) #keep the unrelated pairs without PSE
                unrelated_comp_without_PSE <-c(unrelated_comp_without_PSE, list(paste(paste(ASV2_sampleorigin,ASV2_pair, sep ="_"), 
                                                                                      paste(ASV1_sampleorigin,ASV1_pair, sep ="_") , sep="_"))) # all combinations
              }
            }
          }
        }
      }
    }
  }
  # Create contingency tables 
  cont_table <- data.frame(rbind(c(counts_related_above, counts_related_below),
                                 c(counts_unrelated_above, counts_unrelated_below)))
  colnames(cont_table) <- c("Above_threshold", "Below_threshold")
  rownames(cont_table) <- c("Related_distances", "Unrelated_distances")
  
  return(list("Contingency_table"= cont_table, "Pairs_with_PSEs" = pairs_with_PSE, 
              "Unrelated_pairs_with_PSE" = unrelated_comp_with_PSE, 
              "Unrelated_pairs_without_PSE" = unrelated_comp_without_PSE))
}


# Set working directory
setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/")

# Read metadata and previous results files
metadata <- read.delim("Data_sequencing_method_comparison_n136.txt") 

# 16S V3-V4 results
genera_counts_16S <- read.delim("16S/Mother_Milk-Baby_Feces/Genera_counts.txt", header = T, check.names = F) 
genera_in_pairs_16S <- readRDS("16S/Mother_Milk-Baby_Feces/Genera_counts_pairs.rds")
matching_ASVs_16S <- readRDS("16S/Mother_Milk-Baby_Feces/ASVs_shared_genera.rds")

# 16S-ITS-23S results
genera_counts_16S23S <- read.delim("16S23S/Mother_Milk-Baby_Feces/Genera_counts.txt", header = T, check.names = F) 
genera_in_pairs_16S23S <- readRDS("16S23S/Mother_Milk-Baby_Feces/Genera_counts_pairs.rds")
matching_ASVs_16S23S <- readRDS("16S23S/Mother_Milk-Baby_Feces/ASVs_shared_genera.rds")

#################################################################################
# SNP distances analysis in ASVs of shared genera in maternal-baby oral samples
#################################################################################

# Input: Distance matrices created by **snp-dists** using the MSAs generated in the script MSA_mother_milk_baby_oral.R

#################################################
#A. 16S_23S - Mother Milk vs Baby Feces
#################################################

setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/16S23S/Mother_Milk-Baby_Feces/")

#***************************************
# Load the SNP distance matrix
#***************************************
# List the distance matrices
filenames_16S23S <- list.files(pattern="*.tab", recursive = TRUE)

# Create list of data frame names excluding extension
names_shared_genera <-sub("/.*", "", filenames_16S23S)
names_16S23S <-paste(names_shared_genera , "16S23S", sep = "_")

# Load all files (add 1st column as rownames)
SNP_distances_16S23S <- lapply(filenames_16S23S, function(file) read.delim(file, check.names = FALSE))
names(SNP_distances_16S23S) <- names_16S23S
for (i in 1:length(SNP_distances_16S23S)){
  rownames(SNP_distances_16S23S[[i]]) <- SNP_distances_16S23S[[i]][,1]
  SNP_distances_16S23S[[i]][,1] <- NULL
}

# Get half of the distance matrices (avoid duplicated values)
for (i in 1:length(SNP_distances_16S23S)) {
  SNP_distances_16S23S[[i]][upper.tri(SNP_distances_16S23S[[i]])] <- NA
}


#***************************************************************
# Estimate the PSEs in related vs unrelated mother-baby pairs
#***************************************************************
# Generate all combinations of pair IDs 
mothers <- paste0("Mother_Pair",1:14)
babies <- paste0("Infant_Pair", 1:14)

all_pairs <- expand.grid(mothers, babies)
all_pairs_combined <- paste(all_pairs$Var1, all_pairs$Var2, sep="_")
related_pairs <- paste0( "Mother_Pair", 1:14, "_Infant_Pair",1:14)
unrelated_pairs <- all_pairs_combined[!all_pairs_combined %in% related_pairs]
related_pairs <- sub(".+_", "", related_pairs)

# Generate the empty final data table
final_PSE_data <- data.frame(matrix(nrow=length(all_pairs_combined), ncol=length(names_16S23S)))
colnames(final_PSE_data) <- names_16S23S
rownames(final_PSE_data) <- c(related_pairs, unrelated_pairs)
number_related_pairs <- length(unique(metadata_16S23S$pair_ID))
final_PSE_data$Relatedness <- c(rep("Related", number_related_pairs), 
                                rep("Unrelated", length(all_pairs_combined) - number_related_pairs))

# Number of pairs in which each genus is detected
all_pairs <- table(unlist(genera_in_pairs_16S23S))
total_pairs <- all_pairs[names(all_pairs) %in% names_shared_genera]

# Pairs in which each genus is detected
Staphylococcus_in_pairs <- names(grep("Staphylococcus", genera_in_pairs_16S23S, value = T))
Streptococcus_in_pairs <- names(grep("Streptococcus", genera_in_pairs_16S23S, value = T))


# Get the contigency table and the pairs in which we observe PSEs. Add results to final_PSE_data
PSE_Staphylococcus_16S23S <- exact_PSE_analysis(SNP_distances_16S23S$Staphylococcus_16S23S, matching_ASVs_16S23S)
final_PSE_data[rownames(final_PSE_data) %in% Staphylococcus_in_pairs,"Staphylococcus_16S23S" ] <- "no_PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Staphylococcus_16S23S$Pairs_with_PSEs),"Staphylococcus_16S23S" ] <- "PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Staphylococcus_16S23S$Unrelated_pairs_without_PSE),"Staphylococcus_16S23S" ] <- "no_PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Staphylococcus_16S23S$Unrelated_pairs_with_PSE),"Staphylococcus_16S23S" ] <- "PSE"

PSE_Streptococcus_16S23S <- exact_PSE_analysis(SNP_distances_16S23S$Streptococcus_16S23S, matching_ASVs_16S23S)
final_PSE_data[rownames(final_PSE_data) %in% Streptococcus_in_pairs,"Streptococcus_16S23S" ] <- "no_PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Streptococcus_16S23S$Pairs_with_PSEs),"Streptococcus_16S23S" ] <- "PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Streptococcus_16S23S$Unrelated_pairs_without_PSE),"Streptococcus_16S23S" ] <- "no_PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Streptococcus_16S23S$Unrelated_pairs_with_PSE),"Streptococcus_16S23S" ] <- "PSE"

# Save contingency tables in list
cont_tables_16S23S <- list(PSE_Staphylococcus_16S23S$Contingency_table,
                           PSE_Streptococcus_16S23S$Contingency_table)

names (cont_tables_16S23S) <- names_shared_genera

#################################################
#B. V3V4 - Mother Milk vs Baby Feces 
#################################################

setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/16S/Mother_Milk-Baby_Feces/")

#***************************************
# Load the SNP distance matrix
#***************************************
# List the distance matrices
filenames_V3V4 <- list.files(pattern="*.tab", recursive = TRUE)

# Create list of data frame names excluding extension
names_V3V4 <-sub("/.*", "", filenames_V3V4)
names_V3V4 <-paste(names_V3V4 , "V3V4", sep = "_")

# Load all files (add 1st column as rownames)
SNP_distances_V3V4 <- lapply(filenames_V3V4, function(file) read.delim(file, check.names = FALSE))
names(SNP_distances_V3V4) <- names_V3V4
for (i in 1:length(SNP_distances_V3V4)){
  rownames(SNP_distances_V3V4[[i]]) <- SNP_distances_V3V4[[i]][,1]
  SNP_distances_V3V4[[i]][,1] <- NULL
}

# Get half of the distance matrices (avoid duplicated values)
for (i in 1:length(SNP_distances_V3V4)) {
  SNP_distances_V3V4[[i]][upper.tri(SNP_distances_V3V4[[i]])] <- NA
}

#***************************************************************
# Estimate the PSEs in related vs unrelated mother-baby pairs
#***************************************************************
#Number of pairs in which each genus is detected
all_pairs <- table(unlist(genera_in_pairs_16S))
total_pairs <- all_pairs[names(all_pairs) %in% names_shared_genera]

# Pairs in which each genus is detected
Staphylococcus_in_pairs <- names(grep("Staphylococcus", genera_in_pairs_16S, value = T))
Streptococcus_in_pairs <- names(grep("Streptococcus", genera_in_pairs_16S, value = T))

# Get the contigency table and the pairs in which we observe PSEs. Add results to final_PSE_data
PSE_Staphylococcus_V3V4 <- exact_PSE_analysis(SNP_distances_V3V4$Staphylococcus_V3V4, matching_ASVs_16S)
final_PSE_data[rownames(final_PSE_data) %in% Staphylococcus_in_pairs,"Staphylococcus_V3V4"] <- "no_PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Staphylococcus_V3V4$Pairs_with_PSEs),"Staphylococcus_V3V4"] <- "PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Staphylococcus_V3V4$Unrelated_pairs_without_PSE),"Staphylococcus_V3V4" ] <- "no_PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Staphylococcus_V3V4$Unrelated_pairs_with_PSE),"Staphylococcus_V3V4" ] <- "PSE"

PSE_Streptococcus_V3V4 <- exact_PSE_analysis(SNP_distances_V3V4$Streptococcus_V3V4, matching_ASVs_16S)
final_PSE_data[rownames(final_PSE_data) %in% Streptococcus_in_pairs,"Streptococcus_V3V4" ] <- "no_PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Streptococcus_V3V4$Pairs_with_PSEs),"Streptococcus_V3V4" ] <- "PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Streptococcus_V3V4$Unrelated_pairs_without_PSE),"Streptococcus_V3V4" ] <- "no_PSE"
final_PSE_data[rownames(final_PSE_data) %in% unique(PSE_Streptococcus_V3V4$Unrelated_pairs_with_PSE),"Streptococcus_V3V4" ] <- "PSE"


# Save contingency tables in list
cont_tables_V3V4 <- list(PSE_Staphylococcus_V3V4$Contingency_table,PSE_Streptococcus_V3V4$Contingency_table)
names (cont_tables_V3V4) <- names_shared_genera

# Move relatedness column to the first position
relatedness_column <- final_PSE_data$Relatedness
final_PSE_data <- final_PSE_data[, !(names(final_PSE_data) == "Relatedness")]
final_PSE_data <- data.frame(Relatedness = relatedness_column, final_PSE_data)


########################
#C. Save results 
########################
setwd("~/Desktop/PhD/Projects/Johanne methods paper/Final_paper_version/")

write.table(final_PSE_data,"Result_tables/PSE_table_maternal_milk_infant_feces.txt", sep = "\t", row.names = T, quote = FALSE) 

# 16S V3-V4 results
saveRDS(cont_tables_V3V4, "16S/Mother_Milk-Baby_Feces/Contingency_tables_16S.rds")

# 16S-ITS-23S results
saveRDS(cont_tables_16S23S, "16S23S/Mother_Milk-Baby_Feces/Contingency_tables_16S23S.rds")
