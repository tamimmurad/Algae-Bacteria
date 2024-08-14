#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 01:00:00
#SBATCH -C mem256GB
#SBATCH -J krakenpfpFlow

#Useage: kraken_batch.sh fastaDir krakenOutDir krakenDB condaEnv.
: '
Description:Perform Kraken2 analysis on files located under the directory 
(fastaDir) and output the results to (krakenOutDir).

Kraken 2.1.3 is installed in a local conda environment (meta) which is activated
in the first line. This uses a prebuild database in UPPMAX. 

Input: fasta files directory, a directory for the outputs and location of kraken database (e.g. /sw/data/Kraken2_data/prebuilt/k2_pluspfp_20231009/) and
the name of the conda environment where metaphlan is installed.

Outputs: For each fasta file in the directory, 4 text files are produced: report file showing the # of sequences assigned to each taxonmic ID, 
unclassified sequences file (fasta format), classified sequences (fasta format amending assigned taxonomical ID in the SeqID line) and kraken output file. 
'

fastaDir=$1		#Directory of Fasta files tob analyzed.
krakenOutDir=$2		#Output Directory.
krakenDB=$3		#Kraken database location.
condaEnv=$4		#Conda environment with Kraken installation.

conda activate condaEnv

ls $fastaDir|grep .fasta|while read file
do kraken2  --threads 20  --db $krakenDB --report $krakenOutDir$file.txt --unclassified-out $krakenOutDir$file-uc.txt --classified-out $krakenOutDir$file-c.txt --output $krakenOutDir$file-out.txt  $fastaDir$file
done
