#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 12:00:00
#SBATCH -J metaPairMult

speciesList=$1		#csv file containing species and their reads folder name.
condaEnv=$2		#Output Directory.
prefix=$3
outDir=$4		#Outputs directory

conda activate $condaEnv
awk -F',' '{print $5}' $speciesList|grep $prefix|while read sample
do Out=$outDir$sample/metaphlan/metaphlan_report.txt;
  mkdir $outDir$sample/metaphlan
  In1=$(ls /home/tamimftn/Desktop/algae_bac/data/algae_sequencing_data/raw_data/P26503/$sample/02-FASTQ/220708_A00689_0599_AHH2YTDSX3/*R1_001.fastq.gz)
  In2=$(ls /home/tamimftn/Desktop/algae_bac/data/algae_sequencing_data/raw_data/P26503/$sample/02-FASTQ/220708_A00689_0599_AHH2YTDSX3/*R2_001.fastq.gz) 
  #In1=$(ls  /home/tamimftn/Desktop/algae_bac/data/algae_sequencing_data/raw_data/P28566/ngisthlm00325/files/P28566/$sample/02-FASTQ/230630_A00689_0832_AHTFWGDSX5/*R1_001.fastq.gz);
  #In2=$(ls  /home/tamimftn/Desktop/algae_bac/data/algae_sequencing_data/raw_data/P28566/ngisthlm00325/files/P28566/$sample/02-FASTQ/230630_A00689_0832_AHTFWGDSX5/*R2_001.fastq.gz);

  metaphlan  $In1,$In2 --unclassified_estimation --stat 'avg_g' -t rel_ab_w_read_stats --bowtie2db /home/tamimftn/Desktop/algae_bac/data/databases/metaDB/ --nproc 20 --no_map --input_type fastq -o $Out;
  
  done

