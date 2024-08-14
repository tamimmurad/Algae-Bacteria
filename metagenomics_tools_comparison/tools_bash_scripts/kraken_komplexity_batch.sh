#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 03:00:00
#SBATCH -J kom

#Useage: kraken_komplexity_batch.sh fastaDir outDir  condaEnv.
: '
Description:calculate sequence complexities for each fasta file in fastaDir.

Input: fasta files directory, a directory for the outputs and the name of the conda environment where komplexity tool is installed.

Outputs: For each fasta file in the directory, 1 .tsv file is produced with the complexity of each sequence. 
'


fastaDir=$1
outDir=$2
condaEnv=$3

conda activate $condaEnv
ls $fastaDir |grep '\-c.txt' |while read file
  do cat $fastaDir$file|sed  's/ /$/g'|kz -f  > $outDir$(echo $file|sed 's/.fasta-c.txt//g')-kra-kom.tsv
 done


