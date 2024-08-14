# -*- coding: utf-8 -*-
"""
Spyder Editor
script Name: 
analyse_tools_performance.py

Description:

In this script, the results of all methods of metagenomics analysis developed are used
to generate plots for comparisons in performance in terms of false and true positives.
Additionally, the script provide an insight on the error of the number of assigned reads
to species by the overlap method(best method).  

Procedure:

1- Read the species used tax ids into a dataframe.
2- Provide paths to directories and files of methods' results. Also for the reports of artificial fastas.
3- Define 3 reporting dataframes; one for classified reads report, one for the numbers of fps and one for the number of tps.
4- Populate these dataframe with ready results:Kraken2 and Metaphlan.
5- Iterate over each sample result of Bracken (both thresholded and unthresholded) and get the tps and fps.
6- Iterate over each sample result of the overlap methods and get tps, fps, errors and statistics of reads.
7- Produce plots.

Inputs:

    1- #Paths to artificial fastas and results directories and files.
    2- Species NCBI taxonomy ids table.

          
Outputs:
    
    1- Figure for false positives.
    2- Figure for true positives.
    3- Figure for error in # of reads assigned.


Author: Tamim AlMurad

August 9th 2024.
Version 1.0

"""

import os
import pandas as pd
import matplotlib.pyplot as plt


# =============================================================================
# Files Preperation and Initializations.
# =============================================================================
#Targeted species Tax IDs file.
taxIDsFile = '../species_taxids.tsv'    

#Paths to artificial fastas and results directories and files.
metaResultsFile = '../out_complexity_analysis/metaphlan_results/true_&_false_positives_counts_by_complexity.tsv'    #We can get results of both complexity filtering and no filtering.
krakenResultsFile = '../out_complexity_analysis/kraken_results/true_&_false_positives_counts_by_complexity.tsv'     #Path for Bracken Reports
referencePath = '../out_generated_fasta/'                                                        #Path for simulated fasra reports.
brackenThResultsPath ='../out_tools_results/bracken_100_threshold/'                                             #Path to Bracken with 100 threshold reports.
brackenNoThResultsPath='../out_tools_results/bracken_no_threshold/'                                             #Path to Bracken with no threshold reports.
overlapThResultsPath = '../out_overlap_analysis/overlap_results_100_th/'                                            #Path to overlap with 100 bracken threshold reports.
overlapNoThResultsPath = '../out_overlap_analysis/overlap_results_no_th/'                                           #Path to overlap with no bracken threshold reports.

#Read the targeted tax ids into a dataframe.
taxIDs = pd.read_csv(taxIDsFile,sep='\t')

#Get reference files names (reports genereated by the simulated fasta tool).
refFiles = [report for report in os.listdir(referencePath) if report.endswith('.tsv')]                                                                        #Get the sampling reports file names.

#Read kraken complexity results and seperate FPs and TPs into twp dataframes.
krakenResults = pd.read_csv(krakenResultsFile,sep='\t').set_index('Threshold')
tpKrakenCom=  krakenResults.filter(like='TP',axis=1)
fpKrakenCom=   krakenResults.filter(like='FP',axis=1)                                                          #Get the bracken results file names.

#Read metaphlan complexity results and seperate FPs and TPs into twp dataframes.
metaphlanResults = pd.read_csv(metaResultsFile,sep='\t').set_index('Threshold')                                                                    #Get the Metaphlan results file names.
tpMetaphCom=  metaphlanResults.filter(like='TP',axis=1)
fpMetaphCom=   metaphlanResults.filter(like='FP',axis=1) 

coverages = [str(x) for x in [1,3,5,10,15,20]]                                                               #Coverages' list to be iteratable.
samples =[str(x) for x in [0,1,2,3]]                                                                        #Sampling list to be iteratable-



#A data fame showing total reads and the reads classified by each tool.
classReadsReport =pd.DataFrame()
classReadsReport['coverage'] = coverages
classReadsReport['Ref'] = None
classReadsReport['Bracken'] = None
classReadsReport['Metaphlan'] = None
classReadsReport = classReadsReport.set_index('coverage')

#A data frmae for the number of correctly detected spicies by each tool(True Positives). 
tpReport =pd.DataFrame()
tpReport['coverage'] = coverages
tpReport['Ref'] = [len(taxIDs['taxid'])]*len(coverages)
tpReport['Kraken2'] = list(tpKrakenCom[tpKrakenCom.index==0].loc[0])
tpReport['Kraken2Complex0.8'] = list(tpKrakenCom[tpKrakenCom.index==80].loc[80])
tpReport['BrackenNoTh'] = None
tpReport['BrackenTh100'] = None
tpReport['Metaphlan'] = list(tpMetaphCom[tpMetaphCom.index==0].loc[0])
tpReport['MetaphlanComplex0.7'] = list(tpMetaphCom[tpMetaphCom.index==70].loc[70])
tpReport['OverlapTh100'] = None
tpReport['OverlapNoTh'] = None

tpReport = tpReport.set_index('coverage')

#A data frame for the false positives reported by each tool.
fpReport =pd.DataFrame()
fpReport['coverage'] = coverages
fpReport['Ref'] = 0
fpReport['Kraken2'] = list(fpKrakenCom[fpKrakenCom.index==0].loc[0])
fpReport['Kraken2Complex0.8'] = list(fpKrakenCom[fpKrakenCom.index==80].loc[80])
fpReport['BrackenNoTh'] = None
fpReport['BrackenTh100'] = None
fpReport['Metaphlan'] = list(fpMetaphCom[fpMetaphCom.index==0].loc[0])
fpReport['MetaphlanComplex0.7'] = list(fpMetaphCom[fpMetaphCom.index==70].loc[70])
fpReport['OverlapTh100'] = None
fpReport['OverlapNoTh'] = None

fpReport = fpReport.set_index('coverage')





# =============================================================================
# Fill Bracken Numbers.
# =============================================================================
#Iterate over the coverages.
for cov in coverages:
    #Initalize totals of fps and tps among samples
    tptotalTaxIDTh =0
    tptotalTaxIDNoTh =0
    fptotalTaxIDTh =0
    fptotalTaxIDNoTh =0
    
    #Iterate over the samples.
    for sample in samples:
       
        #Read the bracken result files into dataframes.
        brackenNoTh = pd.read_csv(brackenNoThResultsPath+cov+'-'+sample+'.tsv',sep='\t') 
        brackenTh= pd.read_csv(brackenThResultsPath+cov+'-'+sample+'.tsv',sep='\t')
        
        #Find the number of true positives in each sample.
        sampleTPsTh = len(brackenTh[brackenTh['taxonomy_id'].isin(pd.to_numeric(taxIDs['taxid']))])
        sampleTPsNoTh =len(brackenNoTh[brackenNoTh['taxonomy_id'].isin(pd.to_numeric(taxIDs['taxid']))])
        
        #Sum of samples.
        tptotalTaxIDTh += sampleTPsTh
        tptotalTaxIDNoTh += sampleTPsNoTh
        
        #Find the number of false positives and sum over samples.
        fptotalTaxIDTh += len(brackenTh)-sampleTPsTh
        fptotalTaxIDNoTh += len(brackenNoTh)-sampleTPsNoTh
    
    #Calculate for each coverage the average of fps and tps over samples and assign to the corresponding dataframe.
    tpReport.loc[cov,'BrackenTh100'] = tptotalTaxIDTh/len(samples)      
    tpReport.loc[cov,'BrackenNoTh'] = tptotalTaxIDNoTh/len(samples)
    fpReport.loc[cov,'BrackenTh100'] = fptotalTaxIDTh/len(samples)      
    fpReport.loc[cov,'BrackenNoTh'] = fptotalTaxIDNoTh/len(samples)

# =============================================================================
# Fill Metaphlan and Bracken Overlap Numbers.
# =============================================================================

statIndex = ['mean','var']
tpReadsReportCovTh = pd.DataFrame(columns=coverages,index=statIndex)
tpReadsReportCovNoTh = pd.DataFrame(columns=coverages,index=statIndex)

#Iterate over the coverages.
for cov in coverages:
    tptotalTaxIDTh =0
    tptotalTaxIDNoTh =0
    fptotalTaxIDTh =0
    fptotalTaxIDNoTh =0
    meanTPErrorTh = 0
    meanTPErrorNoTh = 0
    varTPErrorTh = 0
    varTPErrorNoTh = 0
    numOfReads = 0
    brackenClassReads=0
    metaphClassReads=0
   
    #Iterate over the samples.
    for sample in samples:
        
        #Read the reference reports of each sample.
        refReport = pd.read_csv(referencePath+'report'+cov+'-'+sample+'.tsv',sep='\t')
        
        #Corrections to match the species names(Spices should be Spices but it I went through to avoid too many correction)
        
        refReport['Spices']=refReport['Spices'].str.replace('_',' ')
        
        #This piece of code can be commented if the names of species in the simulated fasta reports and tools output are matching.
        refReport.loc[refReport['Spices']=='Chlamydomonas reinhardtii v5','Spices']='Chlamydomonas reinhardtii'
        refReport.loc[refReport['Spices']=='Variovorax sp','Spices']='Variovorax sp. CY25R-8'
        
        #Get the number of reads assigned to each tax ids in the original file.
        tpRefReport = pd.merge(taxIDs, refReport,left_on='species',right_on='Spices')
        neededCol = [x for x in tpRefReport.columns if (x =='taxid' or x=='#of Reads') ]
        tpRefReport = tpRefReport[neededCol]
        
        #Get read stats reports into a dataframe and get totals over the samples.
        statReport = pd.read_csv(overlapNoThResultsPath+cov+'-'+sample+'-reads_stat.tsv',sep='\t')
        numOfReads +=int(statReport.loc[0,'Total_Reads'])
        brackenClassReads +=float(statReport.loc[0,'Bracken%'])
        metaphClassReads +=float(statReport.loc[0,'Metaphlan %'])
        
        #Read the overlap results into dataframes for each sample.
        overlapTh = pd.read_csv(overlapThResultsPath+cov+'-'+sample+'-overlap_result.tsv',sep='\t')
        overlapNoTh = pd.read_csv(overlapNoThResultsPath+cov+'-'+sample+'-overlap_result.tsv',sep='\t')
        
        #Calculate the errors in the number of reads assigned to each taxa for each sample.
        tpReadsReportTh = pd.merge(tpRefReport,overlapTh,left_on='taxid',right_on='taxID')[neededCol+['Avg # of Reads']]
        tpReadsReportTh['Error%'] = 100*abs(tpReadsReportTh['#of Reads']-tpReadsReportTh['Avg # of Reads'])/tpReadsReportTh['#of Reads']
        tpReadsReportNoTh = pd.merge(tpRefReport,overlapTh,left_on='taxid',right_on='taxID')[neededCol+['Avg # of Reads']]
        tpReadsReportNoTh['Error%'] = 100*abs(tpReadsReportTh['#of Reads']-tpReadsReportTh['Avg # of Reads'])/tpReadsReportTh['#of Reads']
        
        #Calculate the mean error and standard deviation among species for each sample and sum over the samples. 
        meanTPErrorTh += tpReadsReportTh['Error%'].mean()
        meanTPErrorNoTh += tpReadsReportTh['Error%'].mean()
        varTPErrorTh += tpReadsReportTh['Error%'].std()
        varTPErrorNoTh+=tpReadsReportNoTh['Error%'].std()
        
        #Find the number of true positives in each sample.
        sampleTPsTh = len(overlapTh[overlapTh['taxonomy_id'].isin(pd.to_numeric(taxIDs['taxid']))])
        sampleTPsNoTh =len(overlapNoTh[overlapNoTh['taxonomy_id'].isin(pd.to_numeric(taxIDs['taxid']))])
        
        #Sum of samples.
        tptotalTaxIDTh += sampleTPsTh
        tptotalTaxIDNoTh += sampleTPsNoTh
        
        #Find the number of false positives and sum over samples.
        fptotalTaxIDTh += len(overlapTh)-sampleTPsTh
        fptotalTaxIDNoTh += len(overlapNoTh)-sampleTPsNoTh
    
    #For each coverage, calculate the average of reads in the original file and the classified over samples and added to the corresponding dataframe.    
    classReadsReport.loc[cov,'Ref']=numOfReads/len(samples)
    classReadsReport.loc[cov,'Bracken']=brackenClassReads/len(samples)
    classReadsReport.loc[cov,'Metaphlan'] = metaphClassReads/len(samples)
    
    #Calculate for each coverage the average error's mean and standard deviation over samples and assign to the corresponding dataframe.
    tpReadsReportCovTh.loc['mean',cov]=meanTPErrorTh/len(samples)
    tpReadsReportCovTh.loc['var',cov]=varTPErrorTh/len(samples)
    tpReadsReportCovNoTh.loc['mean',cov]=meanTPErrorNoTh/len(samples)
    tpReadsReportCovNoTh.loc['var',cov]=varTPErrorNoTh/len(samples)
    
    #Calculate for each coverage the average of fps and tps over samples and assign to the corresponding dataframe.
    tpReport.loc[cov,'OverlapTh100'] = tptotalTaxIDTh/len(samples)      
    tpReport.loc[cov,'OverlapNoTh'] = tptotalTaxIDNoTh/len(samples)
    fpReport.loc[cov,'OverlapTh100'] = fptotalTaxIDTh/len(samples)      
    fpReport.loc[cov,'OverlapNoTh'] = fptotalTaxIDNoTh/len(samples)



# =============================================================================
# Produce Figures.
# =============================================================================

#Define markers and styles.
markers = ['o', 'v', 's', 'x', '^', 'D', 'p', '+', '*']
styles = [f'{marker}-' for marker in markers]

#Plot FPs

plotFP=fpReport.plot(style = styles,use_index=True, ylabel = '# of Species',title = 'False Positives Performance',logy=True,alpha=0.7)
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
plotFP.figure.savefig('False_Positives.svg',bbox_inches='tight')

#Plot TPs
plotTP=tpReport.plot(style = styles,use_index=True, ylabel = '# of Species',title = 'True Positives Performance',alpha=0.7)
plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
plotTP.figure.savefig('True_Positives.svg',bbox_inches='tight',)

#Plot Error stats.
plotStat = plt.figure(figsize=(10,6))
plt.errorbar(tpReadsReportCovTh.columns,tpReadsReportCovTh.loc['mean'],yerr=tpReadsReportCovTh.loc['var'],fmt='--o',color='blue',ecolor='red', capsize=5)
plt.title('Mean and Standard Deviation in Errors of # of Reads among True Positives')
plt.ylabel('Error %')
plt.xlabel('Coverage')
plotStat.figure.savefig('Error_Reads.png')
