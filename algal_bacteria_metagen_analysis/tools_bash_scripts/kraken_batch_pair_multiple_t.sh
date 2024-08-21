#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 04:00:00
#SBATCH -C mem256GB
#SBATCH -J krakenpfpFlow



speciesList=$1		#csv file containing species and their reads folder name.
condaEnv=$2		#Output Directory.
prefix=$3		#Prefix for the parent directory.
outDir=$4		#Outputs directory

conda activate $condaEnv
awk -F',' '{print $5}' $speciesList|grep $prefix| while read sample
   do
   In1=/home/tamimftn/Desktop/algae_bac/data/algae_sequencing_data/raw_data/P26503/$sample/02-FASTQ/220708_A00689_0599_AHH2YTDSX3/*R1_001.fastq.gz
   In2=/home/tamimftn/Desktop/algae_bac/data/algae_sequencing_data/raw_data/P26503/$sample/02-FASTQ/220708_A00689_0599_AHH2YTDSX3/*R2_001.fastq.gz 
   #In1=/home/tamimftn/Desktop/algae_bac/data/algae_sequencing_data/raw_data/P28566/ngisthlm00325/files/P28566/$sample/02-FASTQ/230630_A00689_0832_AHTFWGDSX5/*R1_001.fastq.gz;
   #In2=/home/tamimftn/Desktop/algae_bac/data/algae_sequencing_data/raw_data/P28566/ngisthlm00325/files/P28566/$sample/02-FASTQ/230630_A00689_0832_AHTFWGDSX5/*R2_001.fastq.gz;
   Out=$outDir$sample;
   mkdir $Out
   Out=$Out/kraken/;
   mkdir $Out
   
   echo 'Starting with '$Out
   kraken2 --paired  --threads 20  --db /sw/data/Kraken2_data/prebuilt/k2_pluspfp_20231009/ --report $Out$sample-report.txt --unclassified-out $Out$sample#-uc.txt --classified-out $Out$sample#-c.txt --output $Out$sample-out.txt $In1 $In2;
   done