#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  6 14:30:04 2024
Description:
    This script reads the output relative abundance table produced by the overlap of Metaphlan and Kraken-Bracken
in order to perform analysis and produce figures through defined functions.  

Functions:
    1- produce_stacked_chart(mappingMatrix,aggTaxLevel,algaCellPhen,figType,logScale =None,phylum=None)
    2- produce_pie_chart(mappingMatrix,aggTaxLevel,algaCellPhen,phylum=None)
    3- perform_pca(matrix,aggTaxLevel,uniCellAlg,multiCellAlg)

Inputs:

    1- Relative abundance table as tsv file.
    2- Table of algal sample species and their phenotype (uni/multi).
          
Outputs:

    1- RA (Relative Abundance) stacked chart for all species by bacterial Class (taxonomy level).
    2- RA stacked chart for multicellular algae bacterial species by class.
    3- RA stacked chart for unicellular algae bacterial species by class.
    4- PCA analysis and PCA plot of PC1 and PC2 of RA showing multi/uni cellular algal species.

@author: Tamim AlMurad
"""
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib import cm

RAFile ='/home/tamimftn/Desktop/algae_bac/code/algal_bacteria_metagen_analysis/mapping_results/overlap_RA_Microb_Algae_mapping.tsv'
algSpFile = '/home/tamimftn/Desktop/algae_bac/code/algal_bacteria_metagen_analysis/spieces.csv'
outDir = '/home/tamimftn/Desktop/algae_bac/code/algal_bacteria_metagen_analysis/figures/'

#Load the relative abundance table of all species in the algal samples.
RAinAlgSamps = pd.read_csv(RAFile,sep='\t',index_col=0)

#Remove Eukaryotes species in the results.

bacRAinAlgSamps = RAinAlgSamps[RAinAlgSamps['Kingdom']!='k__Eukaryota']

#Load the algal species with their phelogeny (multi/uni) cellular.
algaCellPhen = pd.read_csv(algSpFile,usecols=['Species','Phenotype'])


#Sorting and clearning.
algaCellPhen.sort_values(by=['Phenotype'],inplace=True)
bacRAinAlgSamps.sort_values(by =[])
algaCellPhen.reset_index(inplace=True)
algaCellPhen.drop(columns=['index'],inplace=True)

#Two seperate lists for multicellular and unicellular algae species  
multiCellAlg=list(algaCellPhen[algaCellPhen['Phenotype']=='Multi']['Species'])
uniCellAlg = list(algaCellPhen[algaCellPhen['Phenotype']=='Uni']['Species'])

#Common bacteria in both algal types
commonBac = bacRAinAlgSamps[(bacRAinAlgSamps['#MultAlgSp'] >0) & (bacRAinAlgSamps['#UniAlgSp'] >0)]

#Bacterial species only present in multicellular algae samples
multiSpecificBac = bacRAinAlgSamps[(bacRAinAlgSamps['#MultAlgSp'] >0) & (bacRAinAlgSamps['#UniAlgSp'] ==0)]

#Bacterial species only present in Unicellular algae samples
uniSpecificBac = bacRAinAlgSamps[(bacRAinAlgSamps['#MultAlgSp'] ==0) & (bacRAinAlgSamps['#UniAlgSp'] >0)]

#Unique bacteria found only in one algal spices.
uniqueBac = bacRAinAlgSamps[(bacRAinAlgSamps['#AlgSpTot'] ==1)]

#Create detection matrix
bacDetInAlgSamps = bacRAinAlgSamps.set_index('Microbial Species')
for bac in bacDetInAlgSamps.index:
    for alg in algaCellPhen['Species']:
        if bacDetInAlgSamps.loc[bac,alg] > 0:
            bacDetInAlgSamps.loc[bac,alg]=1

#Data without Chlamydomonas_callosa.
bacRAinAlgSampsWoCall = bacRAinAlgSamps.drop(columns=['Chlamydomonas_callosa'])
bacDetinAlgSampsWoCall = bacDetInAlgSamps.drop(columns=['Chlamydomonas_callosa'])
algaCellPhenWoCall = algaCellPhen[algaCellPhen['Species']!='Chlamydomonas_callosa']
uniCellAlgWoCall = uniCellAlg.copy()
uniCellAlgWoCall.remove('Chlamydomonas_callosa')


# =============================================================================
#  Function produce_stacked_chart. 
# =============================================================================
'''This function takes as required RA table (mapping matrix of species to samples)
    ,required taxa level to be aggregated, the algal species data frame with 
    their phynotype for this plot, for the figure the algal type. 
    ('Multicellular', 'Unicellular' or 'All'). Additional optional arguments are: 
        1- logScale: boolean (True/False). For the y-axis.
        2- phylum: String ('phylum name'). If the anlysis to include only specific phylum.
    
    The output is a tuble containing:
        1- Legend object for the figure which is also saved as png file in the same directory.
        2- A dataframe of the aggregated RA Table.
'''
def produce_stacked_chart(mappingMatrix,aggTaxLevel,algaCellPhen,figType,outDir,logScale =None,phylum=None):
    
    #If only one bacterial phylum is requested to be checked.
    if phylum:
        #Aggregaring RA of the requested phylum over the specified tax level.
        aggrgated = mappingMatrix[mappingMatrix['Phylum']==phylum].groupby([aggTaxLevel],as_index =False)[list(algaCellPhen['Species'])].sum()
    
    #All bacterial phyla are included.
    else:  
        aggrgated = mappingMatrix.groupby([aggTaxLevel],as_index =False)[list(algaCellPhen['Species'])].sum()
        
    #Sorting by most abundant bacterial species.
    aggrgated['Total']=aggrgated.drop([aggTaxLevel],axis=1).sum(axis=1)
    aggrgated.sort_values(by='Total',inplace=True,ascending=False)
    
    #Transposing so algal species are rows.
    aggrgatedT = aggrgated.drop(['Total'],axis=1).set_index(aggTaxLevel).transpose()
    aggrgatedT =pd.merge(aggrgatedT, algaCellPhen,left_on=aggrgatedT.index,right_on='Species')
    aggrgatedT.set_index('Species',inplace=True)
    
    #Sorting with species names.
    aggrgatedT.sort_index(inplace=True)
    
    #Remove unwanted character from bacterial species taxa.
    aggrgatedT.columns=aggrgatedT.columns.str.replace(r'.__','',regex=True)
    
    #Stacked par plot of bacterial RAs for each specified algal species.
    if figType !='All':
        figure,ax = plt.subplots(figsize=(15,8))
        figure = aggrgatedT.plot.bar(stacked=True,ax=ax,rot = 90,colormap=cm.tab20,logy=logScale,fontsize=20).legend(bbox_to_anchor=(1.0,1.0))
        figure.set_title('Bacterial '+ aggTaxLevel)
        plt.title(figType+' Specific Bacteria Relative Abundance by '+ aggTaxLevel + ' in ' + figType +' Algal Samples',fontsize=20)
        plt.savefig(outDir+figType +'by_algae_bac_by_'+aggTaxLevel+'.png', format='png', dpi=300, bbox_inches='tight')
    
    else:
        
        #Adding empty rows for better plot visualization and seperation of mult and uni cellular.
        finalDF = aggrgatedT[aggrgatedT['Phenotype']=='Multi'].copy()
        emptyDF = pd.DataFrame([[None]*len(finalDF.columns)],columns=finalDF.columns,index=[''])
        finalDF=pd.concat([emptyDF,finalDF])
        finalDF=pd.concat([emptyDF,finalDF])
        finalDF=pd.concat([finalDF,emptyDF])
        finalDF=pd.concat([finalDF,emptyDF])
        finalDF = pd.concat([finalDF,aggrgatedT[aggrgatedT['Phenotype']=='Uni']])
        finalDF=pd.concat([finalDF,emptyDF])
        finalDF=pd.concat([finalDF,emptyDF])
        
        figure,ax = plt.subplots(figsize=(15,8))
        figure = finalDF.plot.bar(stacked=True,ax=ax,rot = 90,colormap=cm.tab20,logy=logScale).legend(bbox_to_anchor=(1.0,1.0))
        figure.set_title('Bacterial '+ aggTaxLevel)
        plt.title('Bacteria Relative Abundance by '+ aggTaxLevel + ' in  '+ figType +' Algal samples',fontsize=20)
        plt.savefig(outDir+figType+'_'+aggTaxLevel+'.png', format='png', dpi=300, bbox_inches='tight')
    
    return figure,aggrgatedT

# =============================================================================
# Function produce_pie_chart
# =============================================================================
def produce_pie_chart(mappingMatrix,aggTaxLevel,algaCellPhen,outDir,phylum=None):
    #If only one phylum is requested to be checked.
    if phylum:
        #Aggregaring RA over the specified tax level.
        aggrgated = mappingMatrix[mappingMatrix['Phylum']==phylum].groupby([aggTaxLevel],as_index =False)[list(algaCellPhen['Species'])].sum()
    #All phyla are included.
    else: 
        #Aggregaring RA over the specified tax level.
        aggrgated = mappingMatrix.groupby([aggTaxLevel],as_index =False)[list(algaCellPhen['Species'])].sum()
    
    #Sorting by most abundant bacterial species.
    aggrgated['Total']=aggrgated.drop([aggTaxLevel],axis=1).sum(axis=1)
    aggrgated.sort_values(by='Total',inplace=True,ascending=False)
        
    #Transposing
    aggrgatedT = aggrgated.set_index(aggTaxLevel).transpose()
    
    aggrgatedT.columns=aggrgatedT.columns.str.replace(r'.__','',regex=True)
    #Totals of bacteria RA in each taxonimc group as per the taxa level in the specified algal species.
    totalBac=aggrgatedT.sum()
    
    #Pi chart showing RA of bacteria over all specified species.
    plt.figure()
    piChart= totalBac.plot.pie(colormap=cm.tab20,labels=['']*len(totalBac),autopct='%1.1f%%')
    plt.legend(totalBac.index,bbox_to_anchor=(1.0,1.0),title='Bacterial '+aggTaxLevel)
    plt.title('Bacteria Overall Relative Abundance by '+ aggTaxLevel + ' in '  +' Algal Samples')
    plt.savefig(outDir+'overall_bac_by_'+aggTaxLevel+'.png', format='png', dpi=300, bbox_inches='tight')
    
    return piChart

# =============================================================================
# Funcrion perform_pca 
# =============================================================================
def perform_pca(matrix,aggTaxLevel,uniCellAlg,multiCellAlg,outDir):
    
    
    #Get relevant data only.
    matrix = matrix[[aggTaxLevel]+uniCellAlg+multiCellAlg].copy()
    #Group by the specified taxa level.
    matrix = matrix.groupby([aggTaxLevel],as_index =False).sum()
    matrix.set_index(aggTaxLevel,inplace=True)
    
    #Transpose.
    matrixT = matrix.transpose()
    #Different colors for uni and multicellular.
    colors = []
    for algSp in matrixT.index:
        if algSp in uniCellAlg:
            colors.append('red')
        elif algSp in multiCellAlg:
            colors.append('blue')
    
    #Perform PCA with and get 2 PCs.
    scalar = StandardScaler()
    matrixTScaled = scalar.fit_transform(matrixT)
    pca=PCA(n_components=2)
    pcs = pca.fit_transform(matrixTScaled)
    
    #Plotting PC1 and PC2 and getting explained variance.
    pcsDf = pd.DataFrame(data=pcs,columns=['PC1','PC2'])
    pcsDf.set_index(matrixT.index,inplace=True)
    plt.figure(figsize=(8,6))
    plt.scatter(pcsDf['PC1'],pcsDf['PC2'],c=colors)
    plt.scatter([], [], c='red', label='Unicellular') 
    plt.scatter([], [], c='blue', label='Multicellular') 
    plt.legend(loc='upper right')
    plt.xlabel('PC1')
    plt.ylabel('PC2')
    plt.title('PCA Analysis of Bacterial Abundance in Uni/Multi Cellular Algae '+ aggTaxLevel + '-wise ')
    plt.savefig(outDir+'PCA_BY_'+aggTaxLevel+'.png', format='png', dpi=300, bbox_inches='tight')
    explained_variance = pca.explained_variance_ratio_
    return explained_variance

# =============================================================================
# Generating Figures.
# =============================================================================
#Bacterial classes RA for every algal sample.
bacRAbyClassFig = produce_stacked_chart(bacRAinAlgSamps, 'Class', algaCellPhen, 'All',outDir=outDir,logScale=False)

#Bacterial classes RA for overall samples.
bacteriaRAPiChart=produce_pie_chart(bacRAinAlgSamps,'Class',algaCellPhen,outDir=outDir,phylum=None)

#Bacterial classes RA for multicellular specific bacteria.
multiSpBacFig = produce_stacked_chart(multiSpecificBac,'Class',algaCellPhen[algaCellPhen['Phenotype']=='Multi'],'Multicellular',outDir=outDir,logScale=True)

#Bacterial classes RA for unicellular specific bacteria.
uniSpecificBacFig = produce_stacked_chart(uniSpecificBac,'Class',algaCellPhen[algaCellPhen['Phenotype']=='Uni'],'Unicellular',outDir=outDir,logScale=True)

#PCA by species for ORA.
varRASpecies = perform_pca(bacRAinAlgSamps, 'Species', uniCellAlg, multiCellAlg,outDir=outDir)

#PCA by Class for ORA.
varRASpecies = perform_pca(bacRAinAlgSamps, 'Class', uniCellAlg, multiCellAlg,outDir=outDir)

#PCA by species for Detection.
varDetSpecies = perform_pca(bacDetInAlgSamps, 'Species', uniCellAlg, multiCellAlg,outDir=outDir+'det_')

#PCA by Class for Detection.
varDetSpecies = perform_pca(bacDetInAlgSamps, 'Class', uniCellAlg, multiCellAlg,outDir=outDir+'det_')

##PCA by species for ORA without Chlamydomonas_Callosa
varRASpeciesWoCall = perform_pca(bacRAinAlgSampsWoCall, 'Species', uniCellAlgWoCall, multiCellAlg,outDir=outDir+'wo_callosa_')

##PCA by species for Det without Chlamydomonas_Callosa
varDetSpeciesWoCall = perform_pca(bacDetinAlgSampsWoCall, 'Species', uniCellAlgWoCall, multiCellAlg,outDir=outDir+'wo_callosa_det')

##PCA by Class for ORA without Chlamydomonas_Callosa
varRASpeciesWoCall = perform_pca(bacRAinAlgSampsWoCall, 'Class', uniCellAlgWoCall, multiCellAlg,outDir=outDir+'wo_callosa_')

totBacSp = len(bacRAinAlgSamps)
mostCommon = bacRAinAlgSamps[bacRAinAlgSamps['#AlgSpTot'] == max(bacRAinAlgSamps['#AlgSpTot'])]['Microbial Species'].values[0]
mostCommonMult = bacRAinAlgSamps[bacRAinAlgSamps['#MultAlgSp'] == max(bacRAinAlgSamps['#MultAlgSp'])]['Microbial Species'].values[0]
mostCommonUni = bacRAinAlgSamps[bacRAinAlgSamps['#UniAlgSp'] == max(bacRAinAlgSamps['#UniAlgSp'])]['Microbial Species'].values[0]
avgNumBac = totBacSp/len(algaCellPhen) 
