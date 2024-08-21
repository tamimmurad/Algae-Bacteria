#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 02:00:00
#SBATCH -J brakenpfpFlow

#Useage: sbatch bracken_batch.sh krakenDir brackenOutDir krakenDB condaEnv.
: '
Description:Perform Bracken analysis on kraken results files located under
 the directory (krakenDir) and output the results to (krakenOutDir). The 
 settings used are 150bps for sequence length, threshold of 100 reads, and specify the results to have only the species level.

 Bracken 2.9 is installed in a conda environment (e.g. meta) which is activated
Input: directory of Kraken results, a directory for the outputs, location of Kraken database and
the name of the conda environment where Bracken is installed.
Outputs: For kraken report file in the directory,  1 report file showing esch detected species and data about the reads.
'
InDir=$1		#Directory of kraken output reports.
condaEnv=$2		#Conda environment with Bracken installation.
krakenDB=$3
conda activate $condaEnv

ls $InDir |while read sample
do 
  out=$InDir$sample/bracken
  mkdir  $out
  bracken -d $krakenDB -i $InDir$sample/kraken/*report.txt -o $out/bracken_report.txt -r 150 -l S -t 100;
done
