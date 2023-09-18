#!/bin/bash

#SBATCH --job-name=fastqc_and_multiqc
#SBATCH --time=0:59:00
#SBATCH --cpus-per-task=6
#SBATCH --output=fastqc_and_multiqc.out
#SBATCH --error=fastqc_and_multiqc.err
#SBATCH --mem=20GB

# This script was written by Johanne Spreckels.


### FastQC

# load fastqc module
echo $PATH
module load FastQC/0.11.5-Java-1.8.0_144

# make directory for fastqc results
mkdir fastqc_results

## run fastqc
# use files from rawdata folder as input and save output in fastqc_results folder
fastqc -t 6 -o ./fastqc_results MY_PATH_TO_16S_SEQUENCING_RAWDATA_FOLDER/*.fastq.gz


### MultiQC

## load multiqc
module load multiqc/1.7-GCCcore-7.3.0-Python-3.7.4-bare

# make directory for multiqc results
mkdir multiqc_results

## run multiqc
# for forward reads only
multiqc fastqc_results/*R1* -o multiqc_results/ -n multiqc_R1

# for reverse reads only
multiqc fastqc_results/*R2* -o multiqc_results/ -n multiqc_R2

