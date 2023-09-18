# ============================================
#  Script to parse quality trimming err files
# ============================================

# This script was written by Ranko Gacesa, Dec 2021

#start interactive session on the cluster:  srun --cpus-per-task=1 --mem=5gb --nodes=1 --qos=interactive --time=01:00:00 --pty bash -i
#load R module:                             ml RPlus
#start R (R version 4.0.3 (2020-10-10):     R

# set wd
# setwd('<enter WD>')

# prepare data frame for output
res <- NULL

# set file name
inF <- "trimmomatic_quality_trimming.err"

# open file for reading
inFConn <- file(inF,open="r")

#  load file lines
inFL <- readLines(inFConn)

# go over it line by line
resOneRow <- NULL # this is one row of output data frame

# initialize
sName <- NULL; rInput <- NULL; rSurvPE <- NULL; rSurvF <- NULL; rSurvR <- NULL; rDropped <- NULL
for (ln in c(1:length(inFL))) {
  l <- inFL[ln]
  # set state based on what is in there
  if (startsWith(x = l,prefix = " -phred33 ")) {
    # dig out the name of file
    lsF <- unlist(strsplit(l,'/'))[10] #note: you might need to change this number to pull out the correct sample name
    lsR <- unlist(strsplit(l,'/'))[19] #note: you might need to change this number to pull out the correct sample name
    sNameF <- gsub('\\.fastq.gz ','',lsF)
    sNameR <- gsub('\\.fastq.gz ','',lsR)
    # add new set of results to end of results table
  } else if (startsWith(x = l,prefix = "Input Read Pairs: ")) {
    ls <- unlist(strsplit(l,' '))
    rInput <- ls[4]  #note: you might need to change this number to pull out the number of input reads
    rSurvPE <- ls[7] #note: you	might need to change this number to pull out the number	of surviving paired reads
    rSurvF <- ls[12] #note: you	might need to change this number to pull out the number	of surviving forward reads
    rSurvR <- ls[17] #note: you	might need to change this number to pull out the number	of surviving reverse reads
    rDropped <- ls[20] #note: you might need to change this number to pull out the number of dropped reads
  } else if (startsWith(x = l,prefix = "TrimmomaticPE: Completed successfully")) {
    oneRow <- data.frame(Sample.Forward=sNameF,
                         Sample.Reverse=sNameR,
                         Input.Read.Pairs=rInput,
                         Survived.Reads.Paired=rSurvPE,
                         Survived.Reads.Forward=rSurvF,
                         Survived.Reads.Reverse=rSurvR,
                         Dropped.Reads=rDropped)
    res <- rbind.data.frame(res,oneRow)
  }
  # 
}


# close file
close(inFConn)

# save output
write.table(res,'trimmomatic_quality_trimming_read_numbers.csv', sep=',', row.names = F)


