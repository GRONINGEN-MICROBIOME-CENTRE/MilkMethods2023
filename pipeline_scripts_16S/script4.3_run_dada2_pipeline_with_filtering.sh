#!/bin/bash

#SBATCH --job-name=dada2_with_filtering
#SBATCH --time=3:59:00
#SBATCH --cpus-per-task=6
#SBATCH --output=dada2_with_filtering.out
#SBATCH --error=dada2_with_filtering.err
#SBATCH --mem=60GB

# This script was written by Johanne Spreckels.


mkdir dada2_output_with_filtering

# load R module
echo $PATH
module load R/3.6.1-foss-2018b-bare

# run adapted dada2 pipeline with filtering out short and long ASVs
Rscript script4.4_dada2_pipeline_with_filtering.R
