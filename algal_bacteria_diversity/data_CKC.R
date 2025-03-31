#****************************************************************************
#R script for identifying bacteria that are associated with multicellularity in algae
#****************************************************************************


#****************************************************************************
#Packages ####
#****************************************************************************
#Load packages
pacman::p_load(phytools,jtools,arm,interactions,ggplot2,tidyr,dplyr,openxlsx,vegan,reshape2,picante,ape,phyloseq,MCMCglmm)

#****************************************************************************
#Data ####
#****************************************************************************
#relative abundance file.
raw_data = read.csv2('data/inputdata/overlap_RA_Microb_Algae_mapping.tsv',sep='\t')


#Algal species information and keep only the uni/multicellular phenotype column.
alg_pheno = read.csv2('data/inputdata/alg_species.csv',sep=',')
alg_pheno = alg_pheno[,c('Species','Phenotype')]
      
#tree files
alg_tree = read.nexus('data/inputdata/mcc.tre')
bac_tree = read.tree('data/inputdata/all_bac_phyliptree.phy')

#****************************************************************************
#Filtering ####
#****************************************************************************
#Remove non-bacterial species in samples and Callosa data.
bac_wide_df = raw_data[raw_data$Kingdom != 'k__Eukaryota',]

#****************************************************************************
#Formatting ####
#****************************************************************************
#* Trees ----
#**** Bacteria ----
#Cleaning bacterial species names.
bac_tree$tip.label=gsub('\'','',bac_tree$tip.label)

#Checking bacterial species names as in the tree.
for (bacSp in bac_tree$tip.label){
  if (!(bacSp %in% bac_wide_df$Microbial.Species)){
    print(bacSp)
  }
}

#Found the following species with different names.
#Phenylobacterium montanum, 2823693
#Tardiphaga alba,340268 
#Sphaerotilus microaerophilus,2914710
#Roseateles sp. XES5,2877940 

#Correcting bacterial species names to match bacterial tree.
bac_wide_df[bac_wide_df$Microbial.TaxID=='2823693',"Microbial.Species"]='Phenylobacterium montanum'
bac_wide_df[bac_wide_df$Microbial.TaxID=='340268',"Microbial.Species"]='Tardiphaga alba'
bac_wide_df[bac_wide_df$Microbial.TaxID=='2914710',"Microbial.Species"]='Sphaerotilus microaerophilus'
bac_wide_df[bac_wide_df$Microbial.TaxID=='2877940',"Microbial.Species"]='Roseateles sp. XES5'

#Format names 
#bac_tree$tip.label = gsub('.','',bac_tree$tip.label)
bac_tree$tip.label = gsub(' ','_',bac_tree$tip.label)

bac_wide_df$Microbial.Species = gsub(' ','_',bac_wide_df$Microbial.Species)

#Remove unwanted information keeping only bac-alg abundance.
bac_wide_df = bac_wide_df[,names(bac_wide_df) %in% c(alg_pheno$Species,'Microbial.Species')]

#Convert abundance to numeric
bac_wide_df = bac_wide_df %>% mutate(across(2:last_col(), ~ as.numeric(.)))

#Check for ultrametricity
plot(bac_tree,cex=0.5) # some oddities with the tree to look into e.g. lots of polytomies and one branch is shorter than the rest
is.binary(bac_tree)

#make ultrametric
cal = makeChronosCalib(bac_tree, node="root", age.max=2000) #Age of root set according to Marin 2016 MBE: https://academic.oup.com/mbe/article/34/2/437/2740734?login=false
bac_tree = chronos(bac_tree, lambda = 10, model = "correlated", calibration = cal, control = chronos.control(dual.iter.max = 100))
class(bac_tree) = "phylo"

#redo node labels
bac_tree = makeNodeLabel(bac_tree)

#Inverse genentic matrices for analyses
inv_bac_tree = inverseA(bac_tree,nodes="TIPS")$Ainv

#**** Algae ----
#Check names
colnames(bac_wide_df)[colnames(bac_wide_df) %in% alg_tree$tip.label == F]

#fix names
alg_tree$tip.label[alg_tree$tip.label == "Lobomonas_francei"] = "Lobomonas_rostrata"

#Vitreochlamys fluviatilis (which was misidentified in the collection, it is probably actually Haematococcus)
#remove
bac_wide_df = bac_wide_df %>% dplyr::select(-Vitreochlamys_fluviatilis)

#Check for ultrametricity
plot(alg_tree,cex=0.5) #non-ultrametric because of rounding errors
alg_tree = force.ultrametric(alg_tree) #fix

#Inverse genentic matrices for analyses
inv_alg_tree = inverseA(alg_tree)$Ainv

#****************************************************************************
#* Association data ----
#****************************************************************************
#create long dataset(each row = pairwise combination of bac and alg): needed for some analyses
long_df = bac_wide_df %>% pivot_longer(cols=Chlamydomonas_asymmetrica:Volvox_tertius,
                                                   values_to ="rel_abund",names_to = "algae") %>%
                          rename(bacteria="Microbial.Species")

#add in phenotype data and a column of presence absence
long_df = long_df %>% mutate(pheno=alg_pheno$Phenotype[match(algae,alg_pheno$Species)],
                             pheno_num = ifelse(pheno=="Multi",1,0),#useful for modelling probability of multicellularity
                             presence = ifelse(rel_abund >0,1,0),
                             bacphy = bacteria, # for modelling phylogeny
                             algphy = algae,
                             algbacphy = paste(algae,bacteria,sep="_")) # for modelling phylogeny

#****************************************************************************
#* Bacteria associations with multicellularity ----
#****************************************************************************
bac_df = long_df %>% group_by(bacteria) %>% filter(rel_abund > 0) %>% 
                                            summarise(n_multi = sum(pheno_num),
                                                      n_algae = sum(presence),
                                                      n_uni = n_algae - n_multi,
                                                      per_multi = n_multi/n_algae)

#* Create phyloseq object ----
rownames(bac_wide_df) = bac_wide_df$Microbial.Species
alg_bac_mat = as.matrix(bac_wide_df %>% dplyr::select(-Microbial.Species))


sam_dat = long_df %>% group_by(algae,pheno) %>% summarise(bacteria_richness = sum(presence)) %>% as.data.frame()
rownames(sam_dat) = sam_dat$algae

algbac_pseq = phyloseq(otu_table(alg_bac_mat,taxa_are_rows = T), phy_tree(bac_tree),sample_data(sam_dat))

#****************************************************************************
#Alpha diversity of bacteria per algae species ####
#****************************************************************************
#**richness ----
alg_df = long_df %>% group_by(algae,pheno) %>% summarise(bacteria_richness = sum(presence))

#**phylogenetic diversity ----
faith  =  picante::pd(samp = t(as.data.frame(otu_table(algbac_pseq))),
                      tree = phy_tree(algbac_pseq), include.root = F) %>% 
  dplyr::select(-SR) %>% 
  mutate(algae=alg_df$algae)

#Add faith's to the other measures
alg_df =  left_join(alg_df,faith)

#****************************************************************************
#Community composition (beta diversity) ####
#****************************************************************************

#******************************
###** Jaccard ----
jaccard  =  as.matrix(phyloseq::distance(algbac_pseq, method = "jaccard", binary = TRUE))



# replace self-comparisons and upper triange with NA
jaccard[upper.tri(jaccard, diag = TRUE)]  =  NA

# melt and filter NA
jaccard_list  =  jaccard %>%  reshape2::melt() %>% filter(!is.na(value))

#**convert jaccard to shared unshared for modelling
jaccard_list <- jaccard_list %>%
  dplyr::select(Var1, Var2, value) %>%
  mutate(shared = NA, unshared = NA, total_bac = NA)

tmp = long_df %>% filter(presence == 1) %>% dplyr::select(bacteria,algae,presence)

# get total number of unique bacteria shared by algae - warning: takes some time
for (i in  1:nrow(jaccard_list)) { 
  t <- tmp %>% 
    filter(algae == jaccard_list$Var1[i] | algae == jaccard_list$Var2[i]) %>% 
    pivot_wider(names_from = "algae", values_from = "presence", values_fill = 0)
  colnames(t) <- c("bacteria", "algae1", "algae2")
  t <- t %>% mutate(shared = if_else(algae1 == algae2, "shared", "Unique"))
  jaccard_list$shared[i] <-  nrow(t %>% filter(shared == "shared"))
  jaccard_list$total_bac[i] <-  nrow(t)
  jaccard_list$unshared[i] <-  jaccard_list$total_bac[i] - jaccard_list$shared[i]
}

jaccard_list = jaccard_list %>% rename(algae1=Var1,
                                       algae2=Var2,
                                       jaccardindex=value)

#******************************
###** Bray curtis ----
bray  =  as.matrix(phyloseq::distance(algbac_pseq, method = "bray"))

# replace self-comparisons and upper triange with NA
bray[upper.tri(bray, diag = TRUE)]  =  NA

# melt and filter NA
bray_list  =  bray %>%  reshape2::melt() %>% filter(!is.na(value)) %>% mutate(ID=paste(Var1,Var2,sep="_"))

#Merge lists Add in phenotype data
alg_betadiv = jaccard_list %>% mutate(ID=paste(algae1,algae2,sep="_"),
                                     bray=bray_list$value[match(ID,bray_list$ID)],
                                     pheno1=alg_pheno$Phenotype[match(algae1,alg_pheno$Species)],
                                     pheno2=alg_pheno$Phenotype[match(algae2,alg_pheno$Species)],
                                     algae1phy=algae1,#terms required for modelling phylogeny
                                     algae2phy=algae2,
                                     pheno = paste(pheno1,pheno2,sep="_"),
                                     pheno = ifelse(pheno == "Uni_Multi","Multi_Uni",pheno))

#**************************************************************************************
#Output R Files ####
#**************************************************************************************
save(algbac_pseq,
     bac_df,
     alg_df,
     long_df,
     alg_bac_mat,
     alg_betadiv,
     alg_tree,
     bac_tree,
     inv_alg_tree,
     inv_bac_tree,
     file ="data/outputdata/bac_alage_data.Rdata")

#**************************************************************************************
#END ####
#**************************************************************************************
