#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 03:00:00
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
krakenDir=$1		#Directory of kraken output reports.
brackenOutDir=$2	#Output Directory.
krakenDB=$3		#Kraken database location.
condaEnv=$4		#Conda environment with Bracken installation.

conda activate $condaEnv

ls $krakenDir |grep -Ev 'uc.txt|c.txt|out.txt'|while read report
do bracken -d $krakenDB -i $krakenDir$report -o $brackenOutDir$report-br.txt -r 150 -l S -t 100;
done
