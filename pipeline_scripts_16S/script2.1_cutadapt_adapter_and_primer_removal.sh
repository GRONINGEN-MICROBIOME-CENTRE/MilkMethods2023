#!/bin/bash

#SBATCH --job-name=cutadapt_adapter_and_primer_removal
#SBATCH --time=02:59:00
#SBATCH --cpus-per-task=6
#SBATCH --output=cutadapt_adapter_and_primer_removal.out
#SBATCH --error=cutadapt_adapter_and_primer_removal.err
#SBATCH --mem=20GB

# This script was written by Johanne Spreckels.


# load cutadapt module
echo $PATH
module load cutadapt/2.6-GCCcore-7.3.0-Python-3.7.4-bare

# make directory for fastq without adapters and without primers
mkdir adapter_removed_fastq
mkdir adapter_and_primer_removed_fastq


### ADAPTER REMOVAL

# adapters used: Nextera Sequencing adapters
# adapter for fw read 5'-3': ACACTCTTTCCCTACACGACGCTCTTCCGATCT
# adapter for rv read 5'-3': GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT

# for each fw and rv read, remove respective adapter from 5' end (-g) and reverse-complement of the other read's adapter from 3' end (-a)

# change to directory with raw fastq files
cd MY_PATH_TO_16S_SEQUENCING_RAWDATA_FOLDER

# run cutadapt to remove adapters from forward reads (R1)
for i in *R1_001.fastq.gz; do cutadapt -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/adapter_removed_fastq/$i $i; done

# run cutadapt to remove adapters from reverse reads (R2)
for i in *R2_001.fastq.gz; do cutadapt -g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTG -o MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/adapter_removed_fastq/$i $i; done


### PRIMER REMOVAL

# primers used:
# primer for fw read 5'-3': CCTACGGGAGGCAGCAG
# primer for rv read 5'-3': GGACTACHVGGGTWTCTAAT

# for each fw and rv read, remove respective primer from 5' end (-g) and reverse-complement of the other read's primer from 3' end (-a)

# change to directory with adapter-removed fastq files
cd MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/adapter_removed_fastq

# run cutadapt to remove primers from forward reads (R1)
for i in *R1_001.fastq.gz; do cutadapt -g CCTACGGGAGGCAGCAG -a ATTAGAWACCCBDGTAGTCC -o MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/adapter_and_primer_removed_fastq/$i $i; done

# run cutadapt to remove primers from reverse reads (R2)
for i in *R2_001.fastq.gz; do cutadapt -g GGACTACHVGGGTWTCTAAT -a CTGCTGCCTCCCGTAGG -o MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/adapter_and_primer_removed_fastq/$i $i; done








