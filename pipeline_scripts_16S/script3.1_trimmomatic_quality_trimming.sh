#!/bin/bash

#SBATCH --job-name=trimmomatic_quality_trimming
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=6
#SBATCH --output=trimmomatic_quality_trimming.out
#SBATCH --error=trimmomatic_quality_trimming.err
#SBATCH --mem=45GB

# This script was written by Johanne Spreckels.


input="MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/adapter_and_primer_removed_fastq"
run_name=$(basename $input)

mkdir -p MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/trimmomatic_quality_trimmed_fastq/$run_name
output="MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/trimmomatic_quality_trimmed_fastq/"


for file in $input/*

do

base=$(basename $file)
name=${base%_L001*}

echo $name >> $run_name.sample.list.trimmomatic_quality_trimming

done

#ml Java
ml Java/1.8.0_144
cat $run_name.sample.list.trimmomatic_quality_trimming | uniq | while read line

do

java -jar /groups/umcg-llnext/tmp01/umcg-jspreckels/tools/Trimmomatic-0.39/trimmomatic-0.39.jar PE \
  -phred33 /$input/$line\_L001_R1_001.fastq.gz /$input/$line\_L001_R2_001.fastq.gz \
   $output/$run_name/$line\_1_paired.fastq.gz $output/$run_name/$line\_1_unpaired.fastq.gz \
   $output/$run_name/$line\_2_paired.fastq.gz $output/$run_name/$line\_2_unpaired.fastq.gz \
   SLIDINGWINDOW:4:25 MINLEN:50

done

rm $run_name.sample.list.trimmomatic_quality_trimming

# move the trimmed files to the correct output folder
# (they are put into a subfolder within the desired output folder)
mv MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/trimmomatic_quality_trimmed_fastq/adapter_and_primer_removed_fastq/* MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/trimmomatic_quality_trimmed_fastq/
# delete the now empty subfolder
rm -r MY_PATH_TO_16S_SEQUENCING_PROCESSING_FOLDER/trimmomatic_quality_trimmed_fastq/adapter_and_primer_removed_fastq/
