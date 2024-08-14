#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 00:30:00
#SBATCH -J overlap

#Useage: sbatch overlap_metaph_bra_batch.sh scriptFile metaphlanDir BrackenDir outDir condaEnv.
: '
Description: Given the results of metaphlan and bracken, get the overlap species in a tsv file along with some other statstics files.
It uses a python script that takes produce the results.

Input: the script name to perform the analysis, metaphlan relative abundance reports directory, Bracken reports directory, output directory 
and the name of the conda environment where python is installed.

Outputs: For each sample, 3 .tsv files:overlap_result file, overlap_stat for the # of species detected by each tool and the 
overlap and reads_stat for the number of reads classified by the tools and the original number of reads-
'

scriptFile=$1
metaphlanDir=$2
brackenDir=$3
outDir=$4
condaEnv=$5

conda activate $condaEnv
ls $metaphlanDir  |while read file
  do python $scriptFile $metaphlanDir$file $brackenDir$file $outDir$(echo $file|sed 's/.tsv/-/g')
 done


