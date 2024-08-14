#!/bin/bash -l
#SBATCH -A naiss2024-5-186
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 00:20:00
#SBATCH -J generateFastas

#A bash script to produce randome reads in fasta files from a set of genomes files. The script uses different coverages and produces 4 files in each coverage. 
coverage=( 1 3 5 10 15 20 )   	#Required differnet coverages
location=$1			#Directory of the genomes fasta files to be sampled.	
readLength=$2			#Read length of the sequences.
fastaOut=$3			#Directory for the output fastafiles.
report=$4			#A report file prefix for each fasta file.

i=0				#Counter
while [ $i -le 3 ]		#To repeat 4 times.
do				
  for cov in ${coverage[@]}
  do 
    python NGS_fasta_generator.py $location $readLength $fastaOut$cov-$i.fasta $report$cov-$i.tsv $cov
  done
  i=$(( $i + 1 ))
done