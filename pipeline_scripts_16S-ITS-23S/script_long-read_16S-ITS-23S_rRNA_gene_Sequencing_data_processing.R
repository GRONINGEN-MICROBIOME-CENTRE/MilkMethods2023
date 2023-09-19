# ================================================
#  DADA2 Classifier for Long Read (Johanne)
# > now actually with 16S, ITS, 23S sequences
# ================================================
# libs
library(dada2)
library(tidyverse)
# Setup
setwd('~/UMCG/DADA2_tests/')
path = '~/UMCG/Johanne_LR_SBiome/sbpipe_2022_10_26_11_38_23_rgacesa_UnknownWindows/3_phase_filter_pals/'
source('helper_functions_16S.R')

#load 16S/23S copy numbers data
GTDBK.16S23S.full.cn <- read.table('../16S_23S_DB/GTDB_Full_207/stats_copynumber.tsv',
                                   sep='\t',header = T)
# 0) Prep data
# ================================
# Consensus FastQC files
fqs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fqs), ".",fixed = T), `[`, 2)

# test plot quality
plotQualityProfile(fqs[[1]])
hist(nchar(getSequences(fqs[1])), 100)

# 0b) remove primers
nop <- file.path(path, "noprimers", basename(fqs))
# > forward primers
F27 <-      "AGRGTTYGATYMTGGCTCAG"
F27Degen <- "AGRRTTYGATYHTDGYTYAG"
# > reverse primers
R1492 <-     "RGYTACCTTGTTACGACTT" # strain-id
R.StrainID <-"AGTACYRHRARGGAANGR"
R.V1V9r   <- "AAGTCGTAACAAGGTADCBSTA"

# modified F27Degen & StrainID [reversed]
# == filtering & primer removal was done by shoreline thingy software ==
# prim <- removePrimers(fqs[1], nop[1], primer.fwd=F27Degen,
#                       primer.rev=R.V1V9r,
#                       orient=TRUE, verbose=TRUE,compress = F,
#                       max.mismatch = 3, allow.indels = F)
# hist(nchar(getSequences(nop[1])), 100)

# 1) Filter and trim
# ================================
# Place filtered files in filtered/ subdirectory
fqs.noprimer <- sort(list.files(path = file.path(path), pattern=".fastq", full.names = TRUE))
sample.names.noprimer <- sapply(strsplit(basename(fqs.noprimer), ".",fixed = T), `[`, 2)

filtfqs <- file.path(path, "filtered", paste0("sfinderOutput.",sample.names.noprimer, ".fastq"))

# standard filter
out <- filterAndTrim(fqs.noprimer, filtfqs,minQ=2,
                     maxN=0, maxEE=5, truncQ=2, rm.phix=F,
                     compress=F, multithread=F,verbose=T,
                     trimLeft = 20,maxLen = 3000,minLen = 1900,
                     trimRight = 0)

head(out)
#hist(nchar(getSequences(filtfqs)), 100)

# 3) Dereplicate
# ================================
filtfqs <- sort(list.files(path = file.path(path, "filtered"), pattern=".fastq", full.names = TRUE))
dereps <- derepFastq(filtfqs, verbose=TRUE)
#   name the derep-class objects by the sample names
#names(dereps) <- sample.names[1]
#length(dereps$uniques)

# Learn Errors
# ================================
err <-learnErrors(dereps, BAND_SIZE=32, multithread=TRUE, 
                  errorEstimationFunction=dada2:::PacBioErrfun,
                  nbases = 100000000)
# > check results
plotErrors(err, nominalQ=TRUE)

# 4) Denoise
# ================================
# pool: F no pooling; T = all pooled; "pseudo" : pseudo pooling
dadas <- dada(dereps, err=err, multithread=TRUE, pool=F,BAND_SIZE=32)

# 6) make ASV sequences:
# ================================
seqtab <- makeSequenceTable(dadas)
# look at table
dim(seqtab)
# inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# 7) remove chimeras
# ================================
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
# compare number of reads of nochimera / chimera
sum(seqtab.nochim)/sum(seqtab)

saveRDS(seqtab.nochim,file='Johanne_16S23S_seqtab_nochim.RDS')

# 7.1) QC
# ================================
# multi sample
qcTable <- cbind(out, sapply(dadas$denoised, getN), rowSums(seqtab.nochim))
# single sample
#qcTable <- cbind(out, getN(dadas), rowSums(seqtab.nochim))
colnames(qcTable) <- c("reads.filtered", "reads.denoised", "reads.nonchimeric")
#rownames(qcTable) <- sample.names
head(qcTable)
write.table(qcTable,'Johanne_16S23S_qcTable.csv',sep=',')

# 8/b) Assign taxonomy (16S-ITS-23S GTDB full DB)
# =============================================================================
taxa.16S23S.v3 <- assignTaxonomy(seqtab.nochim, "../16S_23S_DB/GTDB_Full_207/full_DADA2.fa", 
                              multithread=TRUE,tryRC = T,verbose = T,minBoot = 40)
tst.16S23S.v3 <- as.data.frame(taxa.16S23S.v3); rownames(tst.16S23S.v3) <- NULL; 
head(tst.16S23S.v3,10)
write.table(taxa.16S23S.v3,'Johanne_16S23S_LR_Taxa.csv',sep=',')
write.table(seqtab.nochim,'Johanne_16S23S_LR_ASVs.csv',sep=',')

# get ASVs counts
# =============================================================================
taxa.16S23S.v3.cnt <- getAsvCntTable(taxa = taxa.16S23S.v3,
                                     asvs = seqtab.nochim)
write.table(taxa.16S23S.v3.cnt,'Johanne_16S23S_LR_AsvTaxCounts.csv',sep=',')
# get abundances
taxa.16S23S.v3.ab <- getAbTblFromCntTbl(taxa.16S23S.v3.cnt)
write.table(taxa.16S23S.v3.cnt,'Johanne_16S23S_LR_AsvTaxAbundances.csv',sep=',')

# aggregate to species level
aggRA <- aggregateAsvAbCntTable2(abTbl = taxa.16S23S.v3.ab,level = 'Species')
write.table(aggRA,'Johanne_16S23S_LR_AggregatedTaxAbundances.csv',sep=',')

aggCnts <- aggregateAsvAbCntTable2(abTbl = taxa.16S23S.v3.cnt,level = 'Species')
write.table(aggCnts,'Johanne_16S23S_LR_AggregatedCounts.csv',sep=',')

# ggplot(tstRA,aes(x=Species,y=Rel.Ab,col=Species,fill=Species)) + geom_col() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# add copy numbers
# tstRA <- addCopyNumbers(aggRA,GTDBK.16S23S.full.cn)
# 
# # plot
# ggplot(tstRA,aes(x=Species,y=Rel.Ab.CopyAdj,col=Species,fill=Species)) + geom_col() + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# # save
# write.table(tstRA,'Johanne_LR_mock_v2.csv',sep=',',row.names = F)
