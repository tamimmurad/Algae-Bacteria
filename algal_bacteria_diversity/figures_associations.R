#****************************************************************************
#R script for analysing if bacteria are associated with multicellularity
#****************************************************************************

#****************************************************************************
#Preliminaries ####
#****************************************************************************
#Load packages
pacman::p_load(dplyr,tidyverse,dplyr,ggplot2,cowplot,ggmap,scatterpie,ape,RColorBrewer,ggnewscale,psych,magick,tidytree,ggtree,ggtreeExtra,tidytext,phyloseq,ggmap)

#Load data file.
load("data/outputdata/bac_alage_data.Rdata")
load("./results/model_results.RData")

#Theme settings
t1<-theme(strip.background=element_rect(fill="white"),text = element_text(size=12),title =element_text(size=12),axis.text.y = element_text(size=12,colour="black"),axis.text.x = element_text(size=12,colour="black"),legend.title = element_text(size=12,colour="black"),legend.text = element_text(size = 10),legend.position="top",legend.justification = "top",panel.grid.major = element_blank(), panel.grid.minor=element_blank(),panel.background=element_rect(colour="black",linetype = 1,fill = "white"))

#Colour settings
pheno_col<-c("grey20","aquamarine3")

#****************************************************************************
#1) Do unicellular and multicellular species differ in their bacterial communities ####
#****************************************************************************
dat = alg_df

#Add columns to identify outliers


dat$richness_distance =abs(dat$bacteria_richness - mean(dat$bacteria_richness))
dat$PD_distance =abs(dat$PD - mean(dat$PD))
richness_outliers = dat %>% arrange(desc(richness_distance)) %>% head(6)
diversity_outliers = dat %>% arrange(desc(PD_distance)) %>% head(6)

dat$normPD=0
dat[dat$pheno=='Multi',]$normPD =dat[dat$pheno=='Multi',]$PD/length(alg_df[alg_df$pheno=='Multi',]$pheno)
dat[dat$pheno=='Uni',]$normPD =dat[dat$pheno=='Uni',]$PD/length(alg_df[alg_df$pheno=='Uni',]$pheno)
dat$normPD_distance =abs(dat$normPD - mean(dat$normPD))
normPD_outliers = dat %>% arrange(desc(normPD_distance)) %>% head(6)

#Richness
a <- ggplot(dat, aes(x = pheno, y = bacteria_richness, fill = pheno)) +
  geom_boxplot(position = position_dodge(width = 0.5), outliers = F, alpha = 0.7) +
  geom_jitter(width = 0.1, height = 0, size = 3, shape = 21, alpha = 0.75) +
  geom_text(data = richness_outliers,
            aes(label = algae, y = bacteria_richness),
            position = position_jitter(width = 0.1), vjust = -1, size = 3) +
  ylab("Bacterial Richness \n ") +
  xlab("") +
  labs(fill = "", name = "") + 
  guides(color = "none", fill = "none") +
  scale_color_manual(values = pheno_col) +
  scale_fill_manual(values = pheno_col) +
  t1
a

#Phylogenetic diversity
b=ggplot(dat, aes(x=pheno,y=PD,fill=pheno))+
  geom_boxplot(position = position_dodge(width = 0.5), outliers = F,alpha=0.7)+
  geom_jitter(width = 0.1, height = 0,size=3,shape=21,alpha=0.75) +
  geom_text(data = diversity_outliers,
            aes(label = algae, y = PD),
            position = position_jitter(width = 0.1), vjust = -1, size = 3) +
  ylab("Phylogenetic diversity \n ")+
  xlab("")+
  labs(fill="",name="")+guides(color="none",fill="none")+
  scale_color_manual(values=pheno_col)+
  scale_fill_manual(values=pheno_col)+
  t1
b



ab<-plot_grid(a,b,nrow=1,ncol=2,labels = "AUTO", align = "hv",label_size = 18)

ggsave(ab,file="./figures/fig1.png",width=12,height=10)

#there is an outlier - what is it?
dat %>% filter(pheno == "Uni" & PD>30000) #Haematoccus 

#Figure 1: The number and diversity of bacteria associated with multicellular and unicellular algae ----

#****************************************************************************
set.seed(888)
ord <- ordinate(algbac_pseq, "NMDS", "bray")
 
#Extract ordination scores
library(vegan)
ord_scores <- as.data.frame(vegan::scores(ord, display = "sites"))  # Extract NMDS sample scores
ord_scores$sample_id <- rownames(ord_scores)  # Add sample IDs for labeling

# Calculate the global centroid and distances for all points

outliers <- ord_scores %>%
  mutate(
    centroid_x = mean(NMDS1),  # Global centroid x-coordinate
    centroid_y = mean(NMDS2),  # Global centroid y-coordinate
    distance = sqrt((NMDS1 - centroid_x)^2 + (NMDS2 - centroid_y)^2)  # Euclidean distance
  ) %>%
  arrange(desc(distance)) %>%  # Sort by distance
  slice_head(n = 5)  # Keep the top 3 outliers

dev.off()
# Plot with labeled top 5 outliers
nmds <- plot_ordination(algbac_pseq, ord, color = "pheno", shape = "pheno") +
  geom_point(size = 3, alpha = 0.75) +
  scale_colour_manual(values = pheno_col) +
  stat_ellipse(geom = "polygon", type = "t", alpha = 0.2, level = 0.95, aes(fill = pheno, color = pheno)) +
  geom_text(data = outliers, aes(x = NMDS1, y = NMDS2, label = sample_id), inherit.aes = FALSE,
            hjust = 1, vjust = -0.2, size = 4, color = "black", fontface = "bold") +
  guides(shape = "none") +
  labs(color = "", name = "") + t1
nmds
# Save the plot
ggsave(nmds, file = "./figures/NMDS-bray.png", width = 8, height = 8)

#Your composition plots could go here as part b?

#Figure 2: The composition of bacteria associated with multicellular and unicellular algae ----

#****************************************************************************
#2.a) Are specific bacteria associated with multicellularity ####
#****************************************************************************
#Select bacteria that are found in 5 or more algae species
dat = as.data.frame(bac_df) %>% filter(n_algae > 1) %>% 
                                mutate(multi_cat = ifelse(per_multi == 1.0,"Multi","Both"),
                                                          multi_cat = ifelse(per_multi == 0.00,"Uni",multi_cat)) #to colour code taxa exclusively found in multicellular or unicellular taxa

#Prune tree
missingdata<-bac_tree$tip.label[bac_tree$tip.label %in% dat$bacteria == "FALSE"]
tree<-drop.tip(bac_tree,missingdata)

dat$node<-nodeid(tree,dat$bacteria)
fig_tree <- full_join(tree, dat, by='node')

a<-ggtree(fig_tree, layout="rectangular", size=0.15, open.angle=5)+
  geom_tiplab(size=3,offset=40) +
  geom_tippoint(aes(x=x+20, size=per_multi,fill = multi_cat), shape=21,colour="black") +
  scale_fill_manual(values=c("grey","aquamarine3","black"))+
  guides(shape="none",fill="none") +
  labs(size="% multicellular",name="")+
  xlim(c(0,2000))+
  theme(legend.position=c(.15, .85),legend.text=element_text(size=15),legend.title = element_text(size=15,colour="black"))
a

ggsave(a,file="./figures/fig3.pdf",width=12,height=17)

#Figure 3: The % of multicellular algae species that bacteria are found in. Green dots indicate bacteria that are exclusively found with multicellular algae and black dots are bacteria exlusively found with unicellular taxa. Data was filtered to only include bacteria found with 2 or more algae ----

#****************************************************************************
#2.b) Are there uni/mult specific bacterial families, genera or species?
#****************************************************************************
####################################
#Function Plot_abundance_diff
#####################################
plot_abundance_diff <- function(raw_data, alg_df, taxonomy_level = "Order", figure_name = "fig4.jpg",prefix = 'o__',height=10) {
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  
  # Convert data to numeric
  tax_raw_data <- raw_data %>%
    mutate(across(12:last_col(), ~ as.numeric(as.character(.))))
  
  # Group by taxonomy level and sum across numeric columns
  aggregated_data <- tax_raw_data %>%
    group_by(across(all_of(taxonomy_level))) %>%
    summarise(across(14:last_col(), ~ sum(as.numeric(.x), na.rm = TRUE)))
  
  # Add average ORA and presence % columns
  aggregated_data <- aggregated_data %>%
    rowwise() %>%
    mutate(
      avgMult = mean(c_across(all_of(alg_df$algae[alg_df$pheno == "Multi"])), na.rm = TRUE),
      avgUni = mean(c_across(all_of(alg_df$algae[alg_df$pheno == "Uni"])), na.rm = TRUE),
      countMulti = round(100 * sum(c_across(all_of(alg_df$algae[alg_df$pheno == "Multi"])) > 0, na.rm = TRUE) / length(alg_df[alg_df$pheno == 'Multi', ]$pheno), 1),
      countUni = round(100 * sum(c_across(all_of(alg_df$algae[alg_df$pheno == "Uni"])) > 0, na.rm = TRUE) / length(alg_df[alg_df$pheno == 'Uni', ]$pheno), 1),
      highlight = (countUni == 0 | countMulti == 0)    
    ) %>%
    ungroup() 
  fullAggregated = aggregated_data 
  aggregated_data = aggregated_data %>% rowwise() %>%
    filter((countUni + countMulti) >= 8) #Filter out families with low presence.
  
  # Remove prefix if applicable
  aggregated_data[[taxonomy_level]] <- sub(prefix, "", aggregated_data[[taxonomy_level]])
  fullAggregated[[taxonomy_level]] <- sub(prefix, "", fullAggregated[[taxonomy_level]])
  
  #Taxa to be highlighted
  highlightedTaxa <- aggregated_data %>%
    filter(countUni == 0 | countMulti == 0) %>%
    pull(.data[[taxonomy_level]])  
  
  #Arrange the data in dscending order of presence in unicellular
  aggregated_data <- aggregated_data %>% arrange(countUni) 
  
  #Make the taxonomy level a factor in the dataframe.
  aggregated_data[[taxonomy_level]] <- factor(
    aggregated_data[[taxonomy_level]], 
    levels = aggregated_data[[taxonomy_level]]
  )
  # Reshape data for ORA Plot
  aggregated_data <- aggregated_data %>% mutate(avgUni = -avgUni)
  aggregated_long <- aggregated_data %>%
    pivot_longer(cols = c(avgMult, avgUni), names_to = "Phenotype", values_to = "Abundance") %>%
    mutate(Phenotype = ifelse(Phenotype == "avgMult", "Multi", "Uni"))
  
  
  # Reshape data for presence Plot
  aggregated_data <- aggregated_data %>% mutate(countUni = -countUni)
  aggregated_long_presence <- aggregated_data %>%
    pivot_longer(cols = c(countUni, countMulti), names_to = "Phenotype", values_to = "Presence") %>%
    mutate(Phenotype = ifelse(Phenotype == "countMulti", "Multi", "Uni"))
  
  # Reshape data for presence Plot
  fullAggregated <- fullAggregated %>% mutate(countUni = -countUni)

  
  # Create the plot
  ORA_plot <- ggplot(aggregated_long, aes(y = .data[[taxonomy_level]], x = Abundance, fill = Phenotype)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_manual(values = c("Multi" = "aquamarine3", "Uni" = "black")) +
    labs(
      title = paste("Average Overlap Relative Abundance by Bacterial", taxonomy_level),
      x = "Average Overlap Relative Abundance",
      y = paste("Bacterial", taxonomy_level),
      fill = "Phenotype"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
    
    # Add count annotations for Multi
    geom_text(
      data = aggregated_data,
      aes(
        y = .data[[taxonomy_level]],
        x = avgMult,
        label = paste0(round(avgMult,4), "%")
      ),
      inherit.aes = FALSE,
      color = "darkgreen",
      hjust = -0.2,
      size = 3
    ) +
    
    # Add count annotations for Uni
    geom_text(
      data = aggregated_data,
      aes(
        y = .data[[taxonomy_level]],
        x = avgUni,
        label = paste0(round(-avgUni,4), "%")
      ),
      inherit.aes = FALSE,
      color = "black",
      hjust = 1.2,
      size = 3
    )+
    #Highlight the interesting ones.
    theme(axis.text.y = element_text(face = ifelse(aggregated_data[[taxonomy_level]] %in% highlightedTaxa, "bold", "plain"),
                                     size = ifelse(aggregated_data[[taxonomy_level]] %in% highlightedTaxa, 13, 10)))
  
  ORA_plot
  
  Presence_plot <- ggplot(aggregated_long_presence, aes(y = .data[[taxonomy_level]], x = Presence, fill = Phenotype)) +
    geom_bar(stat = "identity", position = "identity") +
    scale_fill_manual(values = c("Multi" = "aquamarine3", "Uni" = "black")) +
    labs(
      title = paste("Presence in Algal Samples by Bacterial", taxonomy_level),
      x = "% of Samples",
      y = paste("Bacterial", taxonomy_level),
      fill = "Phenotype"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.text.x = element_text(size = 10),
      plot.title = element_text(hjust = 0.5),
      legend.position = "top"
    ) +
    geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
    
    # Add count annotations for Multi
    geom_text(
      data = aggregated_data,
      aes(
        y = .data[[taxonomy_level]],
        x = countMulti,
        label = paste0(countMulti, "%")
      ),
      inherit.aes = FALSE,
      color = "darkgreen",
      hjust = -0.2,
      size = 3
    ) +
    
    # Add count annotations for Uni
    geom_text(
      data = aggregated_data,
      aes(
        y = .data[[taxonomy_level]],
        x = countUni,
        label = paste0(-countUni, "%")
      ),
      inherit.aes = FALSE,
      color = "black",
      hjust = 1.2,
      size = 3
    ) +
    #Highlight the interesting ones.
    theme(axis.text.y = element_text(face = ifelse(aggregated_data[[taxonomy_level]] %in% highlightedTaxa, "bold", "plain"),
                                     size = ifelse(aggregated_data[[taxonomy_level]] %in% highlightedTaxa, 13, 10)))
  Presence_plot
  # Save the plot
  ggsave(Presence_plot, file = paste0("./figures/", paste0(figure_name,'_pre.jpg')), width = 12, height = height)
  ggsave(ORA_plot, file = paste0("./figures/", paste0(figure_name,'_ORA.jpg')), width = 12, height = height)
  
  return(fullAggregated)
}

specificFam = plot_abundance_diff(raw_data,alg_df,taxonomy_level = 'Family',figure_name = 'specific_fa','f__',10)
specificSp  = plot_abundance_diff(raw_data,alg_df,taxonomy_level = 'Species',figure_name = 'specific_sp','s__',20)
SpecificGen = plot_abundance_diff(raw_data,alg_df,taxonomy_level = 'Genus',figure_name = 'specific_ge','g__',15)
specificOrd = plot_abundance_diff(raw_data,alg_df,taxonomy_level = 'Order',figure_name = 'specific_Or','o__',10)
specificCls = plot_abundance_diff(raw_data,alg_df,taxonomy_level = 'Class',figure_name = 'specific_cl','c__',10)
specificPhy = plot_abundance_diff(raw_data,alg_df,taxonomy_level = 'Phylum',figure_name = 'specific_ph','p__',10)
#****************************************************************************
#2.c) Is there statistical significance in bacterial taxa presence in uni vs multi?  
#****************************************************************************
library(stats)

numOfMultAlg = length(alg_df[alg_df$pheno=='Multi',]$pheno) #Get the number of Multicellular samples.
numOfUniAlg = length(alg_df[alg_df$pheno=='Uni',]$pheno) #Get the number of Unicellular samples.

#Function for chisquare 
chisq_taxa = function(preAbs = NULL,totM= numOfMultAlg, totU=numOfUniAlg){

  preAbs$countUni = -preAbs$countUni
  preAbs$Uni = round(preAbs$countUni *  totU/ 100)
  preAbs$Multi = round(preAbs$countMulti *  totM/ 100)
  preAbs$totalUni = totU
  preAbs$totalMulti = totM

  preAbs$countUni = NULL
  preAbs$countMulti = NULL

  p_values <- numeric(nrow(preAbs))
  # Perform Chi-Square test for each bacterial family
  for (i in 1:nrow(preAbs)) {
    # Construct contingency table
    contingency_table <- matrix(
      c(preAbs$Uni[i], preAbs$totalUni[i] - preAbs$Uni[i],
        preAbs$Multi[i], preAbs$totalMulti[i] - preAbs$Multi[i]),
      nrow = 2, byrow = TRUE
    )

    # Chi-square test
    test_result <- chisq.test(contingency_table)

    # Store p-value
    p_values[i] <- test_result$p.value
  }

  # Apply Benjamini-Hochberg (BH) correction
  p_adj <- p.adjust(p_values, method = "BH")

  # Add results back to dataframe
  preAbs$raw_p_value <- p_values
  preAbs$adjusted_p_value <- p_adj
  return(preAbs)
}
#Function for chisquare implementation 2
mult_prop = function(dat=NULL,taxa="taxa",present_1="",total_1="",present_2="",total_2=""){
  dat = as.data.frame(dat) # remove formatting
  results = data.frame(taxa=character(),#setup dataframe to write results to
                       prop_1 = numeric(),
                       prop_2 = numeric(),
                       lcl = numeric(),
                       ucl = numeric(),
                       statistic = numeric(),
                       pvalue = numeric())
  
  for(i in 1:dim(dat)[1]) {
    res = prop.test(x = c(dat[i,present_1], dat[i,present_2]), #run prop test for each taxa
                    n = c(dat[i,total_1], dat[i,total_2]))
    
    res = data.frame(taxa = dat[i,taxa], #extract results
                     prop_1 = res$estimate[1],
                     prop_2 = res$estimate[2],
                     lcl = res$conf.int[1],
                     ucl = res$conf.int[2],
                     statistic = res$statistic,
                     pvalue = res$p.value)
    results = rbind(results,res) #write to dataframe
  }
  
  results$adjusted_pvalue = p.adjust(results$pvalue, method="BH") # calculate adjusted p values using fdr (Benjamini-Hochberg method)
  return(results)
}

#ChiSquare Family wise
#Apply chisquare with imp.1
preAbsF =chisq_taxa(specificFam[, c("Family", "countUni", "countMulti")]) 
#Apply chisquare with imp.2
chiResultsF = mult_prop(dat = preAbsF, taxa = 'Family',present_1 = 'Uni',total_1 = 'totalUni',
                       present_2 = 'Multi',total_2 = 'totalMulti')
#ChiSquare Genus wise
#Apply chisquare with imp.1
preAbsG =chisq_taxa(SpecificGen[, c("Genus", "countUni", "countMulti")]) 
#Apply chisquare with imp.2
chiResultsG = mult_prop(dat = preAbsG, taxa = 'Genus',present_1 = 'Uni',total_1 = 'totalUni',
                         present_2 = 'Multi',total_2 = 'totalMulti')

#ChiSquare Species  wise
#Apply chisquare with imp.1
preAbsS =chisq_taxa(specificSp[, c("Species", "countUni", "countMulti")]) 
#Apply chisquare with imp.2
chiResultsS = mult_prop(dat = preAbsS, taxa = 'Species',present_1 = 'Uni',total_1 = 'totalUni',
                         present_2 = 'Multi',total_2 = 'totalMulti')

#ChiSquare Order  wise
#Apply chisquare with imp.1
preAbsO =chisq_taxa(specificOrd[, c("Order", "countUni", "countMulti")]) 
#Apply chisquare with imp.2
chiResultsO = mult_prop(dat = preAbsO, taxa = 'Order',present_1 = 'Uni',total_1 = 'totalUni',
                        present_2 = 'Multi',total_2 = 'totalMulti')


#Apply chisquare with imp.1
preAbsC =chisq_taxa(specificCls[, c("Class", "countUni", "countMulti")]) 
#Apply chisquare with imp.2
chiResultsC = mult_prop(dat = preAbsC, taxa = 'Class',present_1 = 'Uni',total_1 = 'totalUni',
                        present_2 = 'Multi',total_2 = 'totalMulti')

#Apply chisquare with imp.1
preAbsP =chisq_taxa(specificPhy[, c("Phylum", "countUni", "countMulti")]) 
#Apply chisquare with imp.2
chiResultsP = mult_prop(dat = preAbsP, taxa = 'Phylum',present_1 = 'Uni',total_1 = 'totalUni',
                        present_2 = 'Multi',total_2 = 'totalMulti')

#****************************************************************************
#3) Do bacteria and algae coevolve 1####
#****************************************************************************
dat = long_df %>% filter(presence>0) %>% dplyr::select(bacteria,algae,rel_abund,pheno) %>%
                  mutate(phenoCol = ifelse(pheno == "Multi","aquamarine3","black")) %>% as.data.frame()
bac_tree = collapse.singles(bac_tree)

figBAassoc<-cophylo(bac_tree,alg_tree,rotate=T,assoc=dat)

jpeg(file="figures/coev.jpeg",width = 20, height = 20, units = 'in', res = 300,quality=100)

plot(figBAassoc,link.lwd=dat$rel_abund,link.lty="solid",lwd=2,fsize=0.4)
tiplabels.cophylo(which="right",pch = 21, bg=dat$phenoCol, cex = 1.5)

dev.off()

#****************************************************************************
#4) Do bacteria and algae coevolve? 2####
#****************************************************************************
library(ape)
association_matrix=otu_table(algbac_pseq)
association_matrix <- t(association_matrix)
association_matrix[association_matrix > 0] <- 1 

algal_species_in_matrix <- rownames(association_matrix)
pruned_alg_tree <- keep.tip(alg_tree, algal_species_in_matrix)



algal_dist <- cophenetic(pruned_alg_tree)
bacterial_dist <- cophenetic(bac_tree)




result <- parafit(algal_dist, bacterial_dist, association_matrix, nperm = 999,test.links = TRUE)

algal_species <- rownames(association_matrix)
bacterial_species <- colnames(association_matrix)

link_results <- data.frame(
  Algal_Species = algal_species[result$link.table[,1]],  # Map host index to name
  Bacterial_Species = bacterial_species[result$link.table[,2]],  # Map parasite index to name
  ParaFit_Statistic = result$link.table[,3],  # Statistic values
  P_Value = result$link.table[,4]  # p-values
)


