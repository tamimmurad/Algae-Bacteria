#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 14:41:43 2024


Script Name: 
    analyze_complex_filtering.py

Description: 

Given true species taxonomical ids, the script provides the numbers if true and 
false positives at different sequence complexity thresholds and output the results into a
plot and .tsv file.


Inputs:
    1- A Directory of the komplexity tool results file.
    2- A .tsv file containing the taxonomical ids of the true species.
    3- An output directory to store the results.
    4- A string indicating the results the input are the results of kraken or metaphlan.
          
Outputs:
    1- A plot of the false positives at different thresholds.
    2- A tsv file showing the counts of true and false positives at each threshold and coverage.
    
Special package:
    ete3 : used to get NCBI taxonomy ids.

Author: Tamim AlMurad

August 14th 2024.
Version 1.0
"""




import pandas as pd
import sys
from ete3 import NCBITaxa
import matplotlib.pyplot as plt

#Check if useage is correct.
if len(sys.argv) != 5 or sys.argv[4] not in ['kra','metaph']:
    
    print('\nProvided %d arguments, please review useage.'%len(sys.argv))
    print("Usage: analyze_complex_filtering.py <komplexity_results_directory> <taxids_file>" 
          "<output_directory>  <kra|metaph>")
    exit()

komDir = sys.argv[1]        #Directory of komplexity tool output after classification.
taxIDsFile = sys.argv[2]    #Species taxonomy ids file.
resultsDir= sys.argv[3]     #Output directory for the resutls.
metagenTool = sys.argv[4]   #Can be 'kra' or 'metaph' corresponding to Kraken2 or Metaphlan



#Define thrsholds to test on.
thresholds = list(range(0,100,10))                                                                                 #thresholds of complexity to test at.

#Define coverages used and the samples for each coverage.
coverages = [str(x) for x in [1,3,5,10,15,20]]                                                               #Coverages' list to be iteratable over file names.
samples =[str(x) for x in [0,1,2,3]]                                                                         #Samples list to be iterable over file names.

# =============================================================================
# Get Descendants tax IDs of the species included.
# =============================================================================
taxIDs = pd.read_csv(taxIDsFile,sep='\t')                       # Read the taxids file into a dataframe.
taxIDs['descendants'] = [[]]*len(taxIDs)                        # Define an empty column of lists.        
ncbi = NCBITaxa()                                               # NCBITaxa object from ete3.
#ncbi.update_taxonomy_database                                  # Updating the local NCBI taxonomy database. Commented to avoid that every time.

# Define a function that takes a taxid and returns its descendant taxa
def get_descendants(taxid):                                     
    descendants=ncbi.get_descendant_taxa(taxid)                 # Get descendants taxids.
    descendants.append(taxid)                                   # Add the taxid itself to the descendants.
    return descendants

taxIDs['descendants'] = taxIDs['taxid'].apply(get_descendants)  # Apply the function to each taxid in the 'taxid' column and assign the result to the 'descendants' column 


# =============================================================================
# Analyse TPs and FPs Threshold Coverage-wise.
# =============================================================================
tpCountsPerCov = pd.DataFrame(columns=coverages)                # True positive counts per coverage.
fpCountsPerCov = pd.DataFrame(columns=coverages)                # False positive counts per coverage.

#Iterate over coverages.
for coverage in coverages:
    
    tpCountsPerSample = pd.DataFrame(columns=samples)           # True positive counts in each sample.
    fpCountsPerSample = pd.DataFrame(columns=samples)           # False positive counts in each sample.
    
    #Iterate over samples.
    for sample in samples:
        
        komplexityResults = pd.read_csv(komDir +coverage+'-'+sample+'-'+metagenTool+'-kom.tsv',
                                        sep='\t',names=['seqID','length','kmer','komplexity'])          #Read a complexity results file into a dataframe.
        komplexityResults.sort_values(by=['seqID'],inplace=True)                                        #Sort by sequence id.
        
        komplexityResults[['seqID','taxID']] = komplexityResults['seqID'].str.split('|',expand=True )   #Split the seqID column to get a seperate column for tax IDs.
        
        komplexityResults=komplexityResults[komplexityResults['taxID']!='']                             #Remove the sequences that are not classified to species level (especially for Metaphlan)

        tpThresholds =[]                                                                                #for each complexity threshold lever, stor the number of true positives.
        fpthresholds = []                                                                              #for each complexity threshold lever, stor the number of false positives.
        
        #Iterate over the threshold levels.
        for threshold in thresholds:
            thrKomplexityResults = komplexityResults[komplexityResults['komplexity']>(threshold/100)]  #Remove sequences with complexity under the current threshold value.
            readsPerTaxa = pd.DataFrame(thrKomplexityResults['taxID'].value_counts()).reset_index()     #Get the number of reads for each taxa of the remaining sequences.
            
            readsPerTaxa['species_taxID'] = ''                                                          #A column for species taxid.
            
            #Iterate over the true species taxonomy ids.
            for spTaxid in taxIDs['taxid']:                                                             
                
                print('starting with the run ' + str(coverage) + str(sample) +              # To show progress.
                       ' threshold '+ str(threshold)+' with species '+str(spTaxid))
                
                #Iterate over the detected species taxids and check if it is among the true species.
                for readTaxid in readsPerTaxa['taxID']:
                    if int(readTaxid) in taxIDs[taxIDs['taxid']==spTaxid]['descendants'].iloc[0]:       #Check if a detected taxa is in the true species tax ids or their subspecies.
                        readsPerTaxa.loc[readsPerTaxa[readsPerTaxa['taxID']==readTaxid].index,'species_taxID']=spTaxid      #Add the true species tax id to the species tax id column.
            
            #Count True and False positives and append to the corresponding dataframe.            
            truePositives = readsPerTaxa[readsPerTaxa['species_taxID']!='']                   #Get the true positives.
            tpThresholds.append(len(truePositives['species_taxID'].unique()))                 #Add teh count of true positives detected at the threshold.
            falsePositves =readsPerTaxa[readsPerTaxa['species_taxID']=='']                    #Get all tax IDs that don't match with any of the true species.
            fpSpecies =[]                                                                     #To store only false positive species.
            for taxid in falsePositves['taxID']:                                              #Loop over the false positives taxid.
                
                lineage =  ncbi.get_rank(ncbi.get_lineage(int(taxid)))                        #Get full lineage of the tax ID.
                if 'species' in lineage.values():                                             #Check if the lineage of the tax id is at species level.
                    fpSpecies.append(list(lineage.keys())[list(lineage.values()).index('species')]) #Add the species level tax id in case of a species or deeper tax rank.
                    
                        
            fpthresholds.append(len(set(fpSpecies)))                        #Add the count of unique false positives detected att the threshold.
            
        tpCountsPerSample[sample]=tpThresholds                                                #Add all thresholds true positives counts to the samples DF.   
        fpCountsPerSample[sample]=fpthresholds                                                #Add all thresholds false positives counts to the samples DF.
    
    tpCountsPerSample['average']=tpCountsPerSample.mean(axis=1)                               #Add a column with average true positives in the four samples.
    fpCountsPerSample['average']=fpCountsPerSample.mean(axis=1)                               #Add a column with average false positives in the four samples.          
    
    tpCountsPerCov[coverage]=tpCountsPerSample['average']                                       # The average of samples is the count in each coverage.
    fpCountsPerCov[coverage]=fpCountsPerSample['average']                                       # The average of samples is the count in each coverage.
    
#Produce one data frame for both false and true positives.
tpCountsPerCov.columns=[s + 'TP' for s in coverages]                                            #Add suffix to column names.
fpCountsPerCov.columns=[s + 'FP' for s in coverages]
tpCountsPerCov['Threshold'] = thresholds                                                        #Add thresholds column.
fpCountsPerCov['Threshold'] = thresholds

tpCountsPerCov= tpCountsPerCov.set_index('Threshold')                                           #Set the index to be the threshold levels.
fpCountsPerCov= fpCountsPerCov.set_index('Threshold')
resultstpfp = pd.merge(tpCountsPerCov, fpCountsPerCov,on='Threshold')                           #Merge false and true positives  in one dataframe.
# =============================================================================
# Results: Plot and Table Output
# =============================================================================

resultstpfp.to_csv(resultsDir+'true_&_false_positives_counts_by_complexity.tsv','\t')       # Write the results into a file.

fpCountsPerCov.index=[round(t*0.01,2) for t in thresholds]                                  # Reindexing with complexity threshold out of 1.
fpCountsPerCov.columns=[ column.replace('FP','') for column in fpCountsPerCov.columns]      # Removing FP from columns names.

plt.figure(0)
plot0=fpCountsPerCov.plot(use_index =True, xlabel = 'Seq Complexity threshold',             #Plotting false positives.
                          ylabel = '# of False Positives',style ='*-',
                          title = 'False Positives with Seq Complexity Threshold')
plt.legend(title='Coverage')
plot0.figure.savefig(resultsDir+'fpVsThrPerCov.png')
