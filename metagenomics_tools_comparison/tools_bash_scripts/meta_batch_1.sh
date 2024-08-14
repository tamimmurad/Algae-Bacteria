#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p core 
#SBATCH -n 8
#SBATCH -t 00:02:00
#SBATCH -J metaphlanFlow

#Useage: meta_batch.sh fastaDir metaOutDir metaDB condaEnv.

: '
Description:Perform metaphlan analysis on files located under the directory 
 (fastaDir) and output the results to (metaOutDir).
Metaphlan 4.1 is installed in the conda environment (meta) which is activated
 in the first line.

Input: fasta files directory, a directory for the outputs and location of Metaphlan database (e.g. /home/tamimftn/Desktop/algae_bac/data/databases/metaDB/) and
the name of the conda environment where metaphlan is installed.

Outputs: For each fasta file in the directory, 2 text files are produced: 1 report file showing each classified sequence ID and the taxonomy assigned to it.
In addition a bowtie2out file which can help reduce the time in further analysis. The bowtie files are stored in the same directory of the fasta files. 
'

fastaDir=$1		#Directory of Fasta files tob analyzed.
metaOutDir=$2		#Output Directory.
metaDB=$3		#Metaphlan database location.
condaEnv=$4		#Conda environment with Metaphlan installation.

conda activate $condaEnv
ls $fastaDir |grep .fasta|while read file
do metaphlan $fastaDir$file  --unclassified_estimation --stat 'avg_g' -t reads_map --bowtie2db $metaDB --nproc 8 --input_type fasta -o $metaOutDir$file.txt
   
done
