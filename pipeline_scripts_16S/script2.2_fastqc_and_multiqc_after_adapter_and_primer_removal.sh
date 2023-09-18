#!/bin/bash

#SBATCH --job-name=fastqc_and_multiqc_after_adapter_and_primer_removal
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=6
#SBATCH --output=fastqc_and_multiqc_after_adapter_and_primer_removal.out
#SBATCH --error=fastqc_and_multiqc_after_adapter_and_primer_removal.err
#SBATCH --mem=20GB

# This script was written by Johanne Spreckels.


### FastQC

# load fastqc module
echo $PATH
module load FastQC/0.11.5-Java-1.8.0_144

# make directory for fastqc results
mkdir fastqc_results_after_adapter_and_primer_removal

## run fastqc
# use files from adapter_and_primer_removed_fastq folder as input and save output in fastqc_results_after_adapter_and_primer_removal folder
fastqc -t 6 -o ./fastqc_results_after_adapter_and_primer_removal MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/adapter_and_primer_removed_fastq/*.fastq.gz


### MultiQC

## load multiqc
module load multiqc/1.7-GCCcore-7.3.0-Python-3.7.4-bare

# make directory for multiqc results
mkdir multiqc_results_after_adapter_and_primer_removal

## run multiqc
# for forward reads only
multiqc fastqc_results_after_adapter_and_primer_removal/*R1* -o multiqc_results_after_adapter_and_primer_removal/ -n multiqc_after_adapter_and_primer_removal_R1

# for reverse reads only
multiqc fastqc_results_after_adapter_and_primer_removal/*R2* -o multiqc_results_after_adapter_and_primer_removal/ -n multiqc_after_adapter_and_primer_removal_R2

