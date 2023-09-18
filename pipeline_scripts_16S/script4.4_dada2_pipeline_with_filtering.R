### ========================== DADA2 PIPELINE ========================== ###

## SCRIPT:           DADA2 PIPELINE
## DESCRIPTION:      The DADA2 pipeline creates an amplicon sequence variant (ASV) table from fastq 16S amplicon sequencing files
## AUTHORS:          Alexander Kurilshikov, Johanne Spreckels
## NOTES:            The script is based on the DADA2 pipeline version 1.16, almost default setup, see: https://benjjneb.github.io/dada2/tutorial.html
##                   The SILVA reference database v138.1 is used for taxonomic assignment, see: https://zenodo.org/record/4587955#.YdQRkBPML0o
## DATE OF CREATION: January 2022


### CONTENTS OF THIS FILE

## 0. BEFORE STARTING
## 1. INPUT DETAILS
## 2. FILTERING AND TRIMMING
## 3. LEARN ERROR RATES
## 4. DENOISING / SAMPLE INFERENCE
## 5. MERGE DENOISED PAIRED READS
## 6. CONSTRUCT AMPLICON SEQUENCE VARIANT (ASV) TABLE
## 7. REMOVE CHIMERAS
## 8. OPTIONAL: FILTER SEQUENCES BY LENGTH
## 9. CREATE READ SUMMARY STATISTICS
## 10. ASSIGN TAXONOMY
## 11. SAVE OUTPUT


### ========================== 0. BEFORE STARTING ========================== ###

## Input data requirements
#  The pipeline assumes that the 16S sequencing data meets the following criteria:
#  - Samples are demultiplexed (= 1 fastq file/sample)
#  - Non-biological nucleotides are removed (adapters, linkers, primers etc.) 
#  - For paired-end sequencing data, forward and reverse fastq files contain reads in matched order

## Install dada2
#  - See: https://benjjneb.github.io/dada2/dada-installation.html

## Reference files
#  Download and unzip the latest version of the SILVA reference files (v138.1) from here: https://zenodo.org/record/4587955#.YdQRkBPML0o
#  - silva_nr99_v138.1_train_set.fa.gz
#  - silva_species_assignment_v138.1.fa.gz


### ========================== 1. INPUT DETAILS ========================== ###

print("### ===== 1. STARTING DADA2 PIPELINE ===== ###")

## Provide the pipeline with paths and file names
.libPaths("/groups/umcg-llnext/tmp01/umcg-jspreckels/R_libraries/")
run_dada2_pipeline = function(
  project.folder = "MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/trimmomatic_quality_trimmed_fastq",
  output.folder = "dada2_output_with_filtering/",
  R1.tag = "1_paired.fastq.gz",
  R2.tag = "2_paired.fastq.gz",
  unite.species = "/groups/umcg-llnext/tmp01/umcg-jspreckels/SILVA_reference_database/silva_species_assignment_v138.1.fa",
  unite.ref = "/groups/umcg-llnext/tmp01/umcg-jspreckels/SILVA_reference_database/silva_nr99_v138.1_train_set.fa",
  visualize = T
) {
  
## Load libraries
library(dada2)

## Save list of sequencing files
files = list.files(project.folder)

## Read in the names of the (matched) fastq files
fnFs = paste(sep="/",project.folder,grep(R1.tag,files,value = T))
fnRs = paste(sep="/",project.folder,grep(R2.tag,files,value = T))

if(!all(sub(R1.tag,"",fnFs)==sub(R2.tag,"",fnRs))) {
  stop("Something is wrong with the names of your sequences, please check\n")
} else {
  sample.names = sub(paste0(project.folder,"/"),"",sub(R1.tag,"",fnFs))
}


### ========================== 2. FILTERING AND TRIMMING ========================== ###

print("### ===== 2. FILTERING AND TRIMMING ===== ###")

## Assign filenames for filtered fastq.gz files
filtFs <- file.path(project.folder, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(project.folder, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) = sample.names
names(filtRs) = sample.names

## Filter and trim reads
# Default settings except for truncLen=c(240,160) which is replaced with minLen = 160
# Note: my sequencing data is already quality trimmed before the dada2 pipeline so I only filter out reads <160bp here but I don't trim; --> you can try playing with minLen
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, minLen = 160,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE


### ========================== 3. LEARN ERROR RATES ========================== ###

print("### ===== 3. LEARN ERROR RATES ===== ###")

## Goal: Learn errors from PCR amplification and sequencing, e.g. potential base swaps; the error frequency should decrease with increasing quality

## Train the error models
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


### ========================== 4. DENOISING / SAMPLE INFERENCE ========================== ###

print("### ===== 4. DENOISING - SAMPLE INFERENCE ===== ###")

## Goal: To find unique sequence variants in the data, e.g. remove reads with low abundance and high likelihood of not having the correct sequence

## Run dada2 clustering
# Default settings process each sample independently.
# Pooling information across samples can increase sensitivity to sequence variants that may be present at very low frequencies in multiple samples.
# --> Try dada(..., pool=TRUE) to pool all samples for sample inference.
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)


### ========================== 5. MERGE DENOISED PAIRED READS ========================== ###

print("### ===== 5. MERGE DENOISED PAIRED READS ===== ###")

## Merge paired reads (not overlapping sequences are removed)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)


### ========================== 6. CONSTRUCT AMPLICON SEQUENCE VARIANT (ASV) TABLE ========================== ###

print("### ===== 6. CONSTRUCT AMPLICON SEQUENCE VARIANT TABLE ===== ###")

## Construct ASV table
seqtab <- makeSequenceTable(mergers)

# Structure of the generated ASV table:
# - Columns = Sequences 
# - Rows = Samples
# - Cells contain abundance counts of the sequences in each sample


### ========================== 7. REMOVE CHIMERAS ========================== ###

print("### ===== 7. REMOVE CHIMERAS ===== ###")

## Check for chimeras
# Chimeras are sequences joined together from ≥2 biological sequences, this is an error commonly occuring in amplicon sequencing when closely related sequences are amplified
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)


### ========================== 8. OPTIONAL: FILTER SEQUENCES BY LENGTH ========================== ###

print("### ===== 8. FILTER SEQUENCES BY LENGTH ===== ###")

## Check sequence length distribution
nchars.seq.nochim = nchar(colnames(seqtab.nochim)) #count the number of letters in the colnames (=sequence length)
pdf(paste0(output.folder, "dada2_histogram_seq_length_nochim.pdf"))
hist(nchars.seq.nochim, xlim=c(100,600))
dev.off()
# --> perhaps filter sequences by length, e.g. choose an interval around the median. Current filtering setting is based on what the ASVs where assigned to.

## Filter sequences by length
seqtab.final = seqtab.nochim[,nchars.seq.nochim > 399 & nchars.seq.nochim < 432] #only keep sequences with 400-431bp length

## Show final sequence length distribution
nchars.seq.final = nchar(colnames(seqtab.final)) #count the number of letters in the colnames (=sequence length)
pdf(paste0(output.folder, "dada2_histogram_seq_length_filtered.pdf"))
hist(nchars.seq.final, xlim=c(100,600))
dev.off()


### ========================== 9. CREATE READ SUMMARY STATISTICS ========================== ###

print("### ===== 9. CREATE READ SUMMARY STATISTICS ===== ###")

## Generate read summary statistics
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim), rowSums(seqtab.final))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "final")
rownames(track) <- sample.names

## --> Check read summary statistics, there should be no step in which the majority of reads are lost!


### ========================== 10. ASSIGN TAXONOMY ========================== ###

print("### ===== 10. ASSIGN TAXONOMY ===== ###")

## Assign taxonomy to sequence variants using the SILVA reference files
#table(nchar(getSequences(seqtab.final)))
taxa <- assignTaxonomy(seqtab.final, unite.ref, multithread = TRUE, tryRC = TRUE) #default settings + tryRC added
taxa <- addSpecies(taxa, unite.species) #add species to the taxonomy table

## Structure of the generated dataframe 'taxa':
# - Columns = Taxonomy levels (Kingdom, Phylum, Class, Order, Family, Genus, Species)
# - Rows = Sequence variants


### ========================== 11. SAVE OUTPUT ========================== ###

print("### ===== 11. SAVE OUTPUT ===== ###")

#length.sequences = nchar(rownames(taxa))
write.table(track, file = paste0(output.folder, "dada2_read_statistics_with_filtering.txt"),sep="\t")
write.table(seqtab.final, file = paste0(output.folder, "ASV_table_with_filtering.txt"),sep="\t")
write.table(taxa,file = paste0(output.folder, "ASV_taxa_with_filtering.txt"),sep="\t")

result = list(dada2_read_statistics_with_filtering = track,
              ASV_table_with_filtering = seqtab.final,
              ASV_taxa_with_filtering = taxa)
}
result  = run_dada2_pipeline(visualize = F)
