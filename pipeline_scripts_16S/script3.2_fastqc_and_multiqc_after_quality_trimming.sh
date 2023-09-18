#!/bin/bash

#SBATCH --job-name=fastqc_and_multiqc_after_quality_trimming
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=6
#SBATCH --output=fastqc_and_multiqc_after_quality_trimming.out
#SBATCH --error=fastqc_and_multiqc_after_quality_trimming.err
#SBATCH --mem=20GB

# This script was written by Johanne Spreckels.


### FastQC

# load fastqc module
echo $PATH
module load FastQC/0.11.5-Java-1.8.0_144

# make directory for fastqc results
mkdir fastqc_results_after_trimmomatic_quality_trimming

## run fastqc
# use ‘_paired’ files from trimmomatic_quality_trimmed_fastq/ folder as input
# save output in fastqc_results_after_trimmomatic_quality_trimming/ folder
fastqc -t 6 -o fastqc_results_after_trimmomatic_quality_trimming MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/trimmomatic_quality_trimmed_fastq/*_paired.fastq.gz


### MultiQC

## load multiqc
module load multiqc/1.7-GCCcore-7.3.0-Python-3.7.4-bare

# make directory for multiqc results
mkdir multiqc_results_after_trimmomatic_quality_trimming

## run multiqc
# for forward reads only
multiqc -o multiqc_results_after_trimmomatic_quality_trimming fastqc_results_after_trimmomatic_quality_trimming/*_1_paired* -n multiqc_after_trimmomatic_quality_trimming_R1 --interactive

# for reverse reads only
multiqc -o multiqc_results_after_trimmomatic_quality_trimming fastqc_results_after_trimmomatic_quality_trimming/*_2_paired* -n multiqc_after_trimmomatic_quality_trimming_R2 --interactive


