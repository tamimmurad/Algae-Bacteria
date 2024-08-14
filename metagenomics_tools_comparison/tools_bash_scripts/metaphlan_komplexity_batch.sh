#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J kom

#Useage: metaphlan_komplexity_batch.sh fastaDir outDir condaEnv.
: '
Description:calculate sequence complexities for each classified sequences of metaphlan results.

Input: metaphlan classified sequences fasta files directory, a directory for the outputs and the name of the conda environment where komplexity is installed.

Outputs: For each fasta file in the directory, 1 .tsv file is produced with the complexity of each sequence. 
'


fastaDir=$1
outDir=$2
condaEnv=$3

conda activate $condaEnv
ls $fastaDir |grep 'c.fasta' |while read file
  do cat $fastaDir$file|kz -f  > $outDir$(echo $file|sed 's/-metaph-c.fasta//g')-metaph-kom.tsv
 done


