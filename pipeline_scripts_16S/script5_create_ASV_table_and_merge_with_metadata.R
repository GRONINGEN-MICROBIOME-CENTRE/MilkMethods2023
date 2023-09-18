### ========================== CREATE ASV TABLE AND MERGE WITH METADATA ========================== ###

## SCRIPT:           CREATE ASV TABLE AND MERGE WITH METADATA
## DESCRIPTION:      
## AUTHORS:          Johanne Spreckels
## NOTES:            
## DATE OF CREATION: January 2023

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=02:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.0.3 (2020-10-10):     R

### CONTENTS OF THIS FILE
## 0. IMPORT DATA
## 1. ENSURE CORRECT STRUCTURE OF ROW/COLNAMES AND DATA IN THE DIFFERENT FILES BEFORE MERGING
## 2. MERGE ASV, taxa, AND reads FILES
## 3. SET CHLOROPLASTS AND MITOCHONDRIA TO UNCLASSIFIED
## 4. MERGE d3 WITH METADATA FILE
## 5. SAVE FILES


### ===== 0. IMPORT DATA ===== ###

## import dada2 output files
# ASV table
ASV <- read.table("MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/dada2_output_with_filtering/ASV_table_with_filtering.txt", header=T)
# rownames = sequencing IDs
# colnames = bacterial sequences
# cells contain number of reads matching each sequence

# ASV taxa
taxa <- read.table("MY_PATH_TO_16S_SEQEUNCING_PROCESSING_FOLDER/dada2_output_with_filtering/ASV_taxa_with_filtering.txt", header=T)
# rownames = bacterial sequences
# colnames = Kingdom, Phylum, Class, Order, Family, Genus, Species
# cells contain bug names

# read statistics from quality trimming
reads1 <- read.csv("MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/trimmomatic_quality_trimming_read_numbers.csv", header=T, sep=",")
# rownames = consecutive numbers
# colnames = step in pipeline at which reads were counted
# cells contain number of reads for given sample at that step

# read statistics
reads2 <- read.table("MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/dada2_output_with_filtering/dada2_read_statistics_with_filtering.txt", header=T)
# rownames = sequencing IDs
# colnames = step in pipeline at which reads were counted
# cells contain number of reads for given sample at that step

## import metadata file
meta <- read.csv("MY_PATH_TO_MY_METADATA_FILE/metadata.csv", header=T, sep=",")


### ===== 1. ENSURE CORRECT STRUCTURE OF ROW/COLNAMES AND DATA IN THE DIFFERENT FILES BEFORE MERGING ===== ###

## ensure ASV, reads1, reads2 and meta have matching sequencing IDs as rownames
#ASV
rownames(ASV) <- gsub(".*_22", "22", rownames(ASV)) #adapt this based on the sample names from the sequencer
rownames(ASV) <- gsub("_S.*", "", rownames(ASV))    #adapt this	based on the sample names from the sequencer

#reads1
rownames(reads1) <- reads1$Sample.Forward
rownames(reads1) <- gsub(".*_22", "22", rownames(reads1)) #adapt this	based on the sample names from the sequencer
rownames(reads1) <- gsub("_S.*", "", rownames(reads1))    #adapt this	based on the sample names from the sequencer

#reads2
rownames(reads2) <- gsub(".*_22", "22", rownames(reads2)) #adapt this	based on the sample names from the sequencer
rownames(reads2) <- gsub("_S.*", "", rownames(reads2))    #adapt this	based on the sample names from the sequencer


## change colnames in reads files
colnames(reads1)[3:7] <- c("qc_trimming_n_reads_input_paired", "qc_trimming_n_reads_survived_paired", "qc_trimming_n_reads_survived_fw", "qc_trimming_n_reads_survived_rv", "qc_trimming_n_reads_dropped")
colnames(reads2) <- c("dada2_n_reads_input", "dada2_n_reads_filtered", "dada2_n_reads_denoised_fw", "dada2_n_reads_denoised_rv", "dada2_n_reads_merged", "dada2_n_reads_nochim", "dada2_n_reads_final")


## prepare taxa file for merging with ASV file
# include abbreviation letter for taxonomic level before entries (to make them grep-able by level later)
taxa$Kingdom <- paste0("k__", taxa$Kingdom)
taxa$Phylum <- paste0("p__", taxa$Phylum)
taxa$Class <- paste0("c__", taxa$Class)
taxa$Order <- paste0("o__", taxa$Order)
taxa$Family <- paste0("f__", taxa$Family)
taxa$Genus <- paste0("g__", taxa$Genus)
taxa$Species <- paste0("s__", taxa$Species)

# merge taxonomic levels into 1 label
cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxa$label <- apply(taxa[,cols], 1, paste, collapse="_")

# save file showing all unique taxa detected in the samples
write.table(unique(taxa$label), "MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/dada2_output_with_filtering/detected_taxa.txt", sep="\t", row.names = F)


### ===== 2. MERGE ASV, taxa, AND reads FILES ===== ###

## merge ASV with taxa (by sequencing ids)
d1 <- merge(taxa, as.data.frame(t(ASV)), by="row.names")
d1$taxon_ASV <- paste0(d1$label, "_ASV__", d1$Row.names)
rownames(d1) <- d1$taxon_ASV
d1 <- d1[,10:(ncol(d1)-1)]
# colnames = sequencing ids
# rownames = taxonomic classification + ASV sequence

## merge d1 with reads2 (by sequencing ids)
d2 <- merge(reads2, as.data.frame(t(d1)), by="row.names")
rownames(d2) <- d2$Row.names
d2 <- d2[,-1]
# colnames = read counts (cols 1:7), taxonomic classification + ASV sequence (cols 8:ncol(d2))
# rownames = sequencing ids

## merge d1 with reads1 (by sequencing ids)
d3 <- merge(reads1, d2, by="row.names")
rownames(d3) <- d3$Row.names
d3 <- d3[,-c(1:3,9)]
# colnames = read counts (cols 1:11), taxonomic classification + ASV sequence (cols 12:ncol(d3))
# rownames = sequencing ids

## ensure sample ids are shown in first column
d3$seq_16S_ID <- rownames(d3)
d3 <- d3[,c(20939,1:20938)]  #adapt to place last column first
# colnames = seq_16S_ID (col 1), read counts (cols 2:12), taxonomic classification + ASV sequence (cols 13:ncol(d3))
# rownames = sequencing ids


### ===== 3. SET CHLOROPLASTS AND MITOCHONDRIA TO UNCLASSIFIED ===== ###

## set Chloroplasts and Mitochondria to unclassified
chloro <- grep("Chloroplast", colnames(d3))
mito <- grep("Mitochondria", colnames(d3))
colnames(d3)[c(chloro,mito)] <- gsub(".*_ASV__", "k__NA_p__NA_c__NA_o__NA_f__NA_g__NA_s__NA_ASV__", colnames(d3)[c(chloro,mito)])

## Archaea are kept in the data set
## Unclassified on kingdom and phylum level are kept in the dataset


### ===== 4. MERGE d3 WITH METADATA FILE ===== ###

## merge d3 with lab metadata (by sequencing ids)
d4 <- merge(meta, d3, by="seq_16S_ID")
# colnames = metadata, taxonomic classification + ASV sequence
# rownames = consecutive numbers


### ===== 5. SAVE FILES ===== ###

## save processed 16S data without metadata
write.table(d3, "MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/processed_16S_data.txt", sep="\t", row.names = F)

## save processed 16S data with metadata
write.table(d4, "MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/processed_16S_data_with_metadata.txt", sep="\t", row.names = F)





