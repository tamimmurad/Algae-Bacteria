#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p core 
#SBATCH -n 8
#SBATCH -t 01:20:00
#SBATCH -J metaphlanFlow

#Useage: meta_batch.sh bowtieDir metaRAOutDir metaDB condaEnv.


: '
Description:Perform metaphlan analysis on previously Metaphlan bowtie files and output the results to (metaRAOutDir).
Metaphlan 4.1 is installed in the conda environment (meta) which is activated
 in the first line.

 Input: bowtie directory, a directory for the outputs, location of Metaphlan database (e.g. /home/tamimftn/Desktop/algae_bac/data/databases/metaDB/) and
the name of the conda environment where Metaphlan is installed.

Outputs: For each bowtie file in the directory, 1 tsv file is produced showing the relative abundance of at each taxonomical level of detected species.
'

bowtieDir=$1		#Directory of Fasta files tob analyzed.
metaRAOutDir=$2		#Output Directory.
metaDB=$3		#Metaphlan database location.
condaEnv=$4		#Conda environment with Metaphlan installation.

conda activate $condaEnv
ls $bowtieDir|grep .fasta|while read file
do metaphlan $bowtieDir$file  --unclassified_estimation --stat 'avg_g' -t rel_ab_w_read_stats --bowtie2db $metaDB --nproc 8 --input_type bowtie2out -o $metaOutDir$file.tsv
   
done
