#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 21 14:49:39 2024

Description:
    This script reads the overlap results of Metaphlan and Bracken of the metagenomics samples.
    It produces a matrix that map each microbial species (mainly bacterial) to each algal species sample.
    The mapping is done with the number of reads and the overlap relative abundance (average abundance between the tools).

Inputs:

    1- Algal species data file.
    2- The directory of the overlap results.
    3- A path for the output directory.
          
Outputs:

    1- Mapping table in a .tsv file with read counts.
    2- Mapping table in .tsv file with overlap relative abundance.
  
@author: Tamim AlMurad
"""



import os, sys
import pandas as pd

#Check correct number of argumens.
if len(sys.argv) != 4:
    print("Usage: python map_algae_bacteria.py <overlap_resutlts_dir> <algal_species_file>" 
          "<output_dir>")
    exit()
    

#Assign arguments to variables.
algaeSpeciesFolder = sys.argv[1]    #Directory where the overlap to be found.
algaCellPhenFile = sys.argv[2]      #file with each algal species and its uni/multi cellular phenotype.  
outputDir =sys.argv[3]

algaCellPhen = pd.read_csv(algaCellPhenFile,usecols=['Species','Phenotype'])     #Read a file with each algal species and its uni/multi cellular phenotype.

#Sorting and cleaning.
algaCellPhen.sort_values(by=['Species'],inplace=True)   
algaCellPhen.reset_index(inplace=True)
algaCellPhen.drop(columns=['index'],inplace=True)

#To get the species folders.
algSpecies = os.listdir(algaeSpeciesFolder) 
algSpecies.sort()
#%%
# =============================================================================
# Get all microbes discovered and calculate overlap relatove abundance.
# =============================================================================

print('Starting by reading overlap results \n')     #To show progress
#A dictionary to store the overlap results of each species.
overlapDict = dict.fromkeys(algSpecies,None)

#A matrix with Microbial species as rows and algal species as columns to save # of reads at each.
speciesAlgaeMatrix = pd.DataFrame(columns=['Microbial Species','Microbial TaxID',
                                           'Clade_name','Kingdom','Phylum','Class',
                                           'Order','Family','Genus','Species',
                                           '#MultAlgSp','#UniAlgSp','#AlgSpTot']+algSpecies)

#Lists to gather data from all samples.
micrSpecies = []
micrTaxIDs=[]
cladeName =[]

#Iterate over all algal species samples.
for sp in algSpecies:
    overlapDict[sp] = pd.read_csv(algaeSpeciesFolder+sp+'/'+'overlap/'+'/overlap_result.tsv',sep='\t',index_col=0)
    #Calculate overlap relative abundance for bracken.
    overlapDict[sp]['Norm_Bracken'] = round(overlapDict[sp]['new_est_reads']/sum(overlapDict[sp]['new_est_reads']),6)
    #Calculate overlap relative abundance for Metaphlan.
    overlapDict[sp]['Norm_Metaphlan'] = round(overlapDict[sp]['estimated_number_of_reads_from_the_clade']/sum(overlapDict[sp]['estimated_number_of_reads_from_the_clade']),6)
    #Calculate Average overlap relative abundance over both tools.
    overlapDict[sp]['Norm_Avg'] = round((overlapDict[sp]['Norm_Bracken']+overlapDict[sp]['Norm_Metaphlan'])/2,6)
    #Add sample's detected microbial species to the lists.
    micrSpecies+=list(overlapDict[sp]['name'])
    micrTaxIDs += list(overlapDict[sp]['taxonomy_id'])
    cladeName += list(overlapDict[sp]['#clade_name'])

#Populate matrix with discovered microbial species and drop duplicates.    
speciesAlgaeMatrix['Microbial Species'] = micrSpecies
speciesAlgaeMatrix['Microbial TaxID']=micrTaxIDs
speciesAlgaeMatrix['Clade_name'] =cladeName
speciesAlgaeMatrix = speciesAlgaeMatrix.drop_duplicates('Microbial TaxID')

#Same as above but to store overlap relative abundance.
speciesAlgaeMatrixNorm = speciesAlgaeMatrix.copy()


#%%
# =============================================================================
# Populate matrices with read counts and relative abundances of each microbial species in each sample. 
# =============================================================================

print('Populating Matrices \n')         #To show progress
for bacSp in speciesAlgaeMatrix['Microbial TaxID']:     #Iterate over rows.
    for algSp in list(overlapDict.keys()):              #Iterate over columns.
        
        if bacSp in overlapDict[algSp]['taxonomy_id'].values:   #Check if microbial species was detected in the sample.
            #Assigned correspending values to the matrices.
            speciesAlgaeMatrix.loc[speciesAlgaeMatrix['Microbial TaxID']==bacSp,[algSp]] = overlapDict[algSp].iloc[overlapDict[algSp][overlapDict[algSp]['taxonomy_id']==bacSp].index.values]['Avg # of Reads'].values
            speciesAlgaeMatrixNorm.loc[speciesAlgaeMatrixNorm['Microbial TaxID']==bacSp,[algSp]]=overlapDict[algSp].iloc[overlapDict[algSp][overlapDict[algSp]['taxonomy_id']==bacSp].index.values]['Norm_Avg'].values
        else:
            #Assign 0 when the microbe was not detected in the sample.
            speciesAlgaeMatrix.loc[speciesAlgaeMatrix['Microbial TaxID']==bacSp,[algSp]]=0
            speciesAlgaeMatrixNorm.loc[speciesAlgaeMatrixNorm['Microbial TaxID']==bacSp,[algSp]]=0
#%%

# =============================================================================
# Add to matrices counts in Multi/Uni cellular algal species.
# =============================================================================

print('Adding counts \n')    #To show progress
#Seperate mult and uni cellular algae.
multiCellAlg=algaCellPhen[algaCellPhen['Phenotype']=='Multi']
uniCellAlg = algaCellPhen[algaCellPhen['Phenotype']=='Uni']

#Iterate over all microbial species.
for algSp in speciesAlgaeMatrix['Microbial TaxID']:
    #One by one count the occurence in mutli, uni cellular and total.
    row =speciesAlgaeMatrix[speciesAlgaeMatrix['Microbial TaxID']==algSp].reset_index()
    numUni =0
    numMulti=0
    for sp in multiCellAlg['Species']:
        if row.loc[0,sp] > 0:
            numMulti+=1
    for sp in uniCellAlg['Species']:
        if row.loc[0,sp] > 0:
            numUni+=1
            
    speciesAlgaeMatrix.loc[speciesAlgaeMatrix.loc[speciesAlgaeMatrix['Microbial TaxID']==algSp].index[0],'#MultAlgSp']=numMulti
    speciesAlgaeMatrix.loc[speciesAlgaeMatrix.loc[speciesAlgaeMatrix['Microbial TaxID']==algSp].index[0],'#UniAlgSp']=numUni
    speciesAlgaeMatrixNorm.loc[speciesAlgaeMatrixNorm.loc[speciesAlgaeMatrixNorm['Microbial TaxID']==algSp].index[0],'#MultAlgSp']=numMulti
    speciesAlgaeMatrixNorm.loc[speciesAlgaeMatrixNorm.loc[speciesAlgaeMatrixNorm['Microbial TaxID']==algSp].index[0],'#UniAlgSp']=numUni


speciesAlgaeMatrix['#AlgSpTot']=speciesAlgaeMatrix['#MultAlgSp']+speciesAlgaeMatrix['#UniAlgSp']
speciesAlgaeMatrixNorm['#AlgSpTot']=speciesAlgaeMatrixNorm['#MultAlgSp']+speciesAlgaeMatrixNorm['#UniAlgSp']

#Seperate lineages into columns.
speciesAlgaeMatrix[['Kingdom','Phylum','Class','Order','Family','Genus','Species']] =speciesAlgaeMatrix['Clade_name'].str.split('|',expand=True)
speciesAlgaeMatrixNorm[['Kingdom','Phylum','Class','Order','Family','Genus','Species']] =speciesAlgaeMatrixNorm['Clade_name'].str.split('|',expand=True)

print('Writing results to %s \n'%outputDir)    #To show progress
#Write the results.
speciesAlgaeMatrix.to_csv(outputDir+'Microb_Algae_mapping.tsv','\t')
speciesAlgaeMatrixNorm.to_csv(outputDir + 'overlap_RA_Microb_Algae_mapping.tsv','\t')

