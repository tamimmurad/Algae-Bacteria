# -*- coding: utf-8 -*-
"""
Spyder Editor
script Name: 

Description:

This script reads the results of running Bracken and Metaphlan on the same fasta file, finds the 
overlap (species detected by both tools) and then provides an output of the overlap with some statisitcs.

Inputs:
    1- Metaphlan results file of the relative abundance setting.
    2- Bracken results file at species level.
    3- An output directory to the overlap results.


          
Outputs:
    1- overlap_result.tsv : a file showing the species overlapped and their data.
    2- overlap_stat.tsv: A file showing the number of species detected by each tool and the overlap.
    3- reads_stat.tsv: A file with statistics about the # of reads original, by each tool and the overlap.


Author: Tamim AlMurad

August 14th, 2024
Version 1.0

"""

import sys
import pandas as pd


#Check correct number of argumens.
if len(sys.argv) != 4:
    print("Usage: python analyze_br_meta_overlap.py <metaphlan_file> <bracken_file>" 
          "<overlap_dir>")
    exit()
    
#Assign arguments to variables.
metaResultsFile = sys.argv[1]
brackenResultsFile = sys.argv[2]
overlapDir= sys.argv[3]

# =============================================================================
# Process Bracken Results.
# =============================================================================


brackenResult = pd.read_csv(brackenResultsFile,sep='\t')         #Get bracken result report into a dataframe.

#Stats
brackenClassReads = sum(brackenResult['new_est_reads'])                                               #Total # of classified reads.

#Spices Detection                       
brackenDetectedSp=len(brackenResult['taxonomy_id'])                                          #Get total number of spices detected in this sample.
brackenResult['taxonomy_id']=brackenResult['taxonomy_id'].apply(str)                                  #Convert taxonomy ids to string.

# =============================================================================
# #Process Metaphlan Result.
# =============================================================================

metaResult = pd.read_csv(metaResultsFile,skiprows=6,sep='\t')      #Read the Metaphlan results of this sample into a dataframe.
metaUncReads = metaResult[metaResult['#clade_name']=='UNCLASSIFIED'].loc[0,'estimated_number_of_reads_from_the_clade']

metaResult=metaResult[metaResult['#clade_name'].str.contains('s__')]                            #Get only spices level from the results.
metaResult=metaResult[metaResult['#clade_name'].str.contains('t__')== False]                    #Remove SGB levels from the results.
metaClassReads = sum(metaResult['estimated_number_of_reads_from_the_clade'])                    #The number of classified reads(reads assigned to detected species)
totalReads=  metaClassReads + metaUncReads                                                      #Total reads in the original sample.

metaSpIds=[]                                                                                    #A list for spices level taxonomy id.
for row in metaResult['clade_taxid']:                                                           #Extract the spices level taxonomy ids.
    metaSpIds.append(row.split('|')[len(row.split('|'))-1])                         
metaResult['taxID'] = metaSpIds                                                                 #Add a column containing the spices level taxonomy id to the results dataframe.

    
metaDetectedSp=(len(metaResult['taxID']))                                                 #Get total number of spices detected in this sample.

# =============================================================================
# Get Overlap between Metaphlan and Bracken Results.
# =============================================================================

#Get the overlaped species detected, calculate the average number of reads assigned and count the species.
overlapResults = pd.merge(brackenResult,metaResult, how='inner',left_on='taxonomy_id',right_on='taxID')     
overlapResults['Avg # of Reads'] =(overlapResults['estimated_number_of_reads_from_the_clade']+overlapResults['new_est_reads'])/2
overlapDetSp=(len(overlapResults['taxID']))

#Produce statstics reports for each tool and the overlap
overlapStatDf = pd.DataFrame.from_dict({'Bracken':[brackenDetectedSp],'Metaphlan':[metaDetectedSp],'Overlap':[overlapDetSp]})
readsCountDF = pd.DataFrame.from_dict({'Bracken Classified':[brackenClassReads],'Bracken%':[round(100*brackenClassReads/totalReads,2)],
'Metaphlan Classified':[metaClassReads],'Metaphlan %':[round(100*metaClassReads/totalReads,2)],'Overlap':[sum(overlapResults['Avg # of Reads'])], 'Total_Reads':[totalReads]})

#To show progress.
print('//////////////Printing files for: ' + overlapDir + '////////////')
print('Bracken: '+str(brackenClassReads)+' Meta: '+str(metaClassReads) + ' Total: '+str(totalReads))

#Output
overlapResults.to_csv(overlapDir+'overlap_result.tsv','\t')
overlapStatDf.to_csv(overlapDir+'overlap_stat.tsv','\t')
readsCountDF.to_csv(overlapDir+'reads_stat.tsv','\t')

