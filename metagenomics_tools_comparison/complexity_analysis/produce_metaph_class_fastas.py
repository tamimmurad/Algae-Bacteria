#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 10:12:35 2024

Script Name: 
    produce_metaph_class_fastas.py

Description: 

Given a file with metaphlan output seqs mapping, extract the corresponding sequences from the fasta file
and produce a new fastafile of only classified sequences. This file will have the original sequence id suffixed with the 
assigned species taxid. The script works on multiple files that are stored in two directory:
    1- Directory for the metaphlan sequence mapping results of the fasta files.
    2- Directory of the fasta files themselves.
Resulting fasta files are stored in the same directory of metaphlan results.
Procedure:



Inputs:

    1- A path to the metaphlan sequence mapping results files.
    2- A path to the fasta files directory.

          
Outputs:
    
    1- For each fasta file, a new fasta file containing only the sequences classified by Metaphlan. 

Author: Tamim AlMurad

August 9th 2024.
Version 1.0
"""

import pandas as pd, os,sys
from Bio import SeqIO

#Check if useage is correct.
if len(sys.argv) != 4:
    
    print('\nProvided %d arguments, please review useage.'%len(sys.argv))
    print("Usage: produce_metaph_class_fastas.py <fastas_directory> <metaphlan_directory>" 
          "<output_directory> ")
    exit()

#Get input and output  directories.
fastaFilesDir = sys.argv[1]
metaphlanSeqDir = sys.argv[2]
outputDit = sys.argv[3]


#To make sure we only get fasta files and not report files.
fastaFiles = [fasta for fasta in os.listdir(fastaFilesDir) if fasta.endswith('.fasta')]  
                                                                  #Get the sampling reports file names.


#Iterate over fastafiles.
for file in fastaFiles:
        
        #Read the sequences into a dictionary and then convert it to a dataframe.
        fastaDic ={rec.id : str(rec.seq) for rec in SeqIO.parse(fastaFilesDir+file,'fasta')}
        fastaDF=pd.DataFrame.from_dict(list(fastaDic.items()))
        fastaDF.columns=['seqID','Seq']
        
        print(file+' is parsed')    #For progress and debugging.
        
        #Read the corresponding metaphlan result file and remove what metaphlan added for sequence ids.
        metaResults = pd.read_csv(metaphlanSeqDir+file+'.txt',sep='\t',skiprows=4) #The first 4 rows contains information that are not needed.
        metaResults['#read_id']=metaResults['#read_id'].str.split('__',expand=True)[0]
        metaResults=metaResults.rename(columns={'#read_id':'seqID'})
        metaResults['taxID'] = metaResults['NCBI_taxlineage_ids'].str.split('|',expand =True).iloc[:,6] #To get only species taxonomy id.
     

        #Get the classified sequences ids and actual sequence and suffix sequences ids with the taxonomy ids of corresponding species. 
        metaResults=pd.merge(metaResults,fastaDF,on='seqID')
        metaResults['seqID'] =metaResults['seqID']+'$metaphlan:taxid|'+metaResults['taxID']

        #Write the results into a new fasta file.
        with open(outputDit+file.replace('.fasta','-metaph-c.fasta'),'w') as handle:
            for index, seq in metaResults.iterrows():
                handle.write('>'+seq['seqID']+'\n'+seq['Seq']+'\n')
  