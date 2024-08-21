#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:10:00
#SBATCH -J overlap

#Useage: sbatch overlap_metaph_bra_batch.sh scriptFile InDir condaEnv outDir 
: '
Description: Given the results of metaphlan and bracken, get the overlap species in a tsv file along with some other statstics files.
It uses a python script that takes produce the results.

Input: the script name to perform the analysis, directory where metaphlan and bracken results, output directory 
and the name of the conda environment where python is installed.

Outputs: For each sample, 3 .tsv files:overlap_result file, overlap_stat for the # of species detected by each tool and the 
overlap and reads_stat for the number of reads classified by the tools and the original number of reads-
'

scriptFile=$1
InDir=$2
condaEnv=$3
i=1
conda activate $condaEnv
ls $InDir  |while read species
  do 
    echo 'Starting with species number' $i
    resultDir=$InDir$species/overlap
    mkdir $resultDir
    python $scriptFile $InDir$species/metaphlan/metaphlan_report.txt $InDir$species/bracken/bracken_report.txt $resultDir/
    i=$((i+1))
 done


