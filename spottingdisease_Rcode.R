# ============================================================
'R code for Microbiome analysis

Christina Pavloudi
cpavloudi@gwu.edu
https://cpavloud.github.io/mysite/

	Copyright (C) 2023 Christina Pavloudi
  
    This script is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
  
    This script is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.'

# =============================================================


############################LOAD LIBRARIES #######################################

# List of packages needed
.packages = c("vegan", "ecodist", "GGally", "BiocManager", "ggplot2", "tidyverse", "RColorBrewer", "DECIPHER", "microViz", "ape")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

# Load packages into session 
lapply(.packages, require, character.only=TRUE)

BiocManager::install("phyloseq")
install.packages("phyloseq")
library(phyloseq)

# Define a default theme for ggplot graphics
theme_set(theme_bw()) 

########################################################################################
############################# 16S rRNA amplicons #######################################
########################################################################################

############################ Preparation of Data #######################################
#import the OTU table (or else biotic data)
bac <- read.csv("seq_table_noblank.csv", sep = ",", header=TRUE, row.names = 1)

#import the taxonomy table
taxonomybac <- read.csv("seq_Taxonomy_silva_merged.csv", sep = ",", header=FALSE, row.names = 1, na.strings=c("","NA"))

colnames(taxonomybac) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#check where there are NA values in the taxonomy data frame
colSums(is.na(taxonomybac))

taxonomy <- taxonomybac

#fill in Phylum column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Phylum[i])){
      taxonomy$Phylum[i]=taxonomy$Kingdom[i]
    }
  if (sum(is.na(taxonomy$Phylum))==0) {
    break
  }
}

#fill in Class column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Class[i])){
      taxonomy$Class[i]=taxonomy$Phylum[i]
    }
  if (sum(is.na(taxonomy$Class))==0) {
    break
  }
}

#fill in Order column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Order[i])){
      taxonomy$Order[i]=taxonomy$Class[i]
    }
  if (sum(is.na(taxonomy$Order))==0) {
    break
  }
}

#fill in Family column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Family[i])){
      taxonomy$Family[i]=taxonomy$Order[i]
    }
  if (sum(is.na(taxonomy$Family))==0) {
    break
  }
}

#fill in Genus column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Genus[i])){
      taxonomy$Genus[i]=taxonomy$Family[i]
    }
  if (sum(is.na(taxonomy$Genus))==0) {
    break
  }
}

#fill in Species column
repeat {
  for (i in 1:nrow(taxonomy))
    if (is.na(taxonomy$Species[i])){
      taxonomy$Species[i]=taxonomy$Genus[i]
    }
  if (sum(is.na(taxonomy$Species))==0) {
    break
  }
}


taxonomybac <- taxonomy

#check where there are NA values in the taxonomy data frame and in the biotic data
colSums(is.na(taxonomybac))
colSums(is.na(bac))

#convert the taxonomy data from data frame to matrix
taxonomy_matrix_bac <- as.matrix(taxonomybac)

# prepare the object for the phyloseq object
TAX_BAC = tax_table(taxonomy_matrix_bac)

#convert the biotic data from data frame to matrix
biotic_matrix_bac <- as.matrix(bac)

#tranpose biotic data for the calculation of diversity indices
biotic_presabs <- decostand(bac, method = "pa")
biotic_presabs_matrix <- as.matrix(biotic_presabs)
OTU_PR = otu_table(biotic_presabs_matrix, taxa_are_rows = TRUE)

#prepare the object for the phyloseq object
OTU_BAC = otu_table(biotic_matrix_bac, taxa_are_rows = TRUE)

#import the metadata of the samples
metadata_physeq <- read.csv("metadata_edited.csv", header=TRUE, row.names = 1) 
# prepare the objects for the phyloseq object
META = sample_data(metadata_physeq)
head(META)

#load the tree file
TREE <- read_tree("seqs_aligned.fa.treefile")
#rename tree tip labels 
TREE$tip.label <- taxa_names(TAX_BAC)

#check if the names are correct and the same
setequal(taxa_names(TAX_BAC), taxa_names(TREE))

######################## PHYLOSEQ analysis #######################################

# combine them all to create the phyloseq object
physeq_bac = phyloseq(OTU_BAC, TAX_BAC, META, TREE)
physeq_PR = phyloseq(OTU_PR, TAX_BAC, META, TREE)
physeq = physeq_bac

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq) == 0)
sum(taxa_sums(physeq) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
# Remove samples that are now empty, ie. that have no counts
physeq <- prune_samples(sample_sums(physeq) > 0, physeq)

#Check if there are ASVs with no counts and how many there are
any(taxa_sums(physeq_PR) == 0)
sum(taxa_sums(physeq_PR) == 0)
#Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
physeq_PR <- prune_taxa(taxa_sums(physeq_PR) > 0, physeq_PR)
# Remove samples that are now empty, ie. that have no counts
physeq_PR <- prune_samples(sample_sums(physeq_PR) > 0, physeq_PR)

#get the data frame from the phyloseq object
pd <- psmelt(physeq)
pd_PR <- psmelt(physeq_PR)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd[,c("Phylum")])))
HowManyPhylaPR <- length(levels(as.factor(pd_PR$Phylum)))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)
PhylaPalettePR = getPalette(HowManyPhylaPR)

# Changing order of levels
pd$Health_status <- factor(pd$Health_status,      # Reordering group factor levels
                           levels = c("Healthy", "Diseased"))
sample_data(physeq)$Health_status <- factor(sample_data(physeq)$Health_status, 
                                            levels = c("Healthy", "Diseased"), ordered = TRUE)
levels(sample_data(physeq)$Health_status)
pd_PR$Health_status <- factor(pd_PR$Health_status,      # Reordering group factor levels
                              levels = c("Healthy", "Diseased"))


#merge the OTUs at the Phylum level
physeq_merged_Phylum <- tax_glom(physeq, "Phylum")
ps0 <- transform_sample_counts(physeq_merged_Phylum, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "16S_OTU_TAX_Merged_Phylum.csv")

#merge the OTUs at the Class level
physeq_merged_Class <- tax_glom(physeq, "Class")
ps0 <- transform_sample_counts(physeq_merged_Class, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "16S_OTU_TAX_Merged_Class.csv")

#merge the OTUs at the Species level
physeq_merged_Species <- tax_glom(physeq, "Species")
ps0 <- transform_sample_counts(physeq_merged_Species, function(x) x)
# Extract abundance matrix from the phyloseq object
OTU_merged = as(otu_table(ps0), "matrix")
# Coerce to data.frame
OTU_merged_df = as.data.frame(OTU_merged)
OTU_merged_df <- tibble::rownames_to_column(OTU_merged_df, "OTU")
# Extract taxonomy matrix from the phyloseq object
TAX_merged = as(tax_table(ps0), "matrix")
# Coerce to data.frame
TAX_merged_df = as.data.frame(TAX_merged)
TAX_merged_df <- tibble::rownames_to_column(TAX_merged_df, "OTU")
#Merge OTU and taxonomy data frames
OTU_TAX_merged <- merge(OTU_merged_df,TAX_merged_df,by = "OTU")
#save the final phyto otu table with the taxonomies
write.csv(OTU_TAX_merged, "16S_OTU_TAX_Merged_Species.csv")



######################### LefSe etc. ###########################################
# 
# # combine them all to create the phyloseq object
# physeq_bac = phyloseq(OTU_BAC, TAX_BAC, META)
# physeq = physeq_bac
# 
# #Check if there are ASVs with no counts and how many there are
# any(taxa_sums(physeq) == 0)
# sum(taxa_sums(physeq) == 0)
# #Remove ASVs/OTUs/taxa that are empty, ie. that have no counts
# physeq <- prune_taxa(taxa_sums(physeq) > 0, physeq)
# # Remove samples that are now empty, ie. that have no counts
# physeq <- prune_samples(sample_sums(physeq) > 0, physeq)

#For the LefSe
library(microbiomeMarker); packageVersion("microbiomeMarker"); citation("microbiomeMarker")
library(MicrobiotaProcess); packageVersion("MicrobiotaProcess"); citation("MicrobiotaProcess")
library(patchwork)
# for the kruskal_test and wilcox_test
library(coin)

physeq_C1 <- subset_samples(physeq, C1.Diseased_microbiome_vs_Healthy_microbiome=="Yes")
physeq_C6 <- subset_samples(physeq, C6.Lesion_tissue_vs_Lesion_swabs=="Yes")
physeq_C8 <- subset_samples(physeq, C8.Sea_urchin_tissue._body_wall_._CF=="Yes")

#LefSe on ASV level
OTU_C1_ASV <- run_lefse(
  physeq_C1,
  group = "Health_status",
  #subgroup = "Sample_type",
  taxa_rank = "none", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 3.0,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)


#LefSe on ASV level
OTU_C1_Species <- run_lefse(
  physeq_C1,
  group = "Health_status",
  #subgroup = "Sample_type",
  taxa_rank = "Species", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)


#LefSe on ASV level
OTU_C6_ASV <- run_lefse(
  physeq_C6,
  group = "Category",
  #subgroup = "Sample_type",
  taxa_rank = "none", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 3.7,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)


#LefSe on ASV level
OTU_C6_Species <- run_lefse(
  physeq_C6,
  group = "Category",
  #subgroup = "Sample_type",
  taxa_rank = "Species", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)


#LefSe on ASV level
OTU_C8_ASV <- run_lefse(
  physeq_C8,
  group = "Group",
  subgroup = "Sample_type",
  taxa_rank = "none", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 3.5,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)


#LefSe on ASV level
OTU_C8_Species <- run_lefse(
  physeq_C8,
  group = "Health_status",
  subgroup = "Sample_type",
  taxa_rank = "Species", #"none" means perform differential analysis on the original taxa e.g., OTU or ASV).
  transform = c("identity"), #return the original data without any transformation (default).
  norm = "CPM", # do not normalize.
  kw_cutoff = 0.05,
  lda_cutoff = 2,
  bootstrap_n = 30,
  bootstrap_fraction = 2/3,
  wilcoxon_cutoff = 0.05,
  multigrp_strat = FALSE,
  strict = "0",
  sample_min = 7,
  only_same_subgrp = FALSE,
  curv = FALSE
)


write.csv(marker_table(OTU_C1_ASV), "lefse_OTU_C1_ASV.csv")
write.csv(marker_table(OTU_C1_Species), "lefse_OTU_C1_Species.csv")
write.csv(marker_table(OTU_C6_ASV), "lefse_OTU_C6_ASV.csv")
write.csv(marker_table(OTU_C6_Species), "lefse_OTU_C6_Species.csv")
write.csv(marker_table(OTU_C8_ASV), "lefse_OTU_C8_ASV.csv")
write.csv(marker_table(OTU_C8_Species), "lefse_OTU_C8_Species.csv")


################# Detailed graphs for each comparison ##########################

####################### C1 ####################################################

physeq <- physeq_C1

#get the data frame from the phyloseq object
pd <- psmelt(physeq)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)

# Changing order of levels
pd$Health_status <- factor(pd$Health_status,      # Reordering group factor levels
                           levels = c("Healthy", "Diseased"))
sample_data(physeq)$Health_status <- factor(sample_data(physeq)$Health_status, 
                                            levels = c("Healthy", "Diseased"), ordered = TRUE)
levels(sample_data(physeq)$Health_status)


alpha_meas = c("Observed")
p_alpha_meas <- plot_richness(physeq, "Health_status", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas+
  geom_boxplot(data=p_alpha_meas$data, aes(x=Health_status, y=value, fill=Health_status), fatten=NULL)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), 
               width= 0.75, linetype= "solid")+
  scale_fill_manual(values= c("#0077b6", "#ae2012"))+
  geom_jitter(width=0, size= 3.5, color= "black")+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  labs(title= "Phylogenetic Diversity", y= "Diversity Index", x="Group")+
  scale_x_discrete(limits= c("Diseased", "Healthy"))

res.aov <- aov(value ~ Health_status, data= p_alpha_meas$data)  
summary(res.aov)

#Alpha diversity surface microbiomes HvD 
alpha_meas = c("Chao1")
p_alpha_meas <- plot_richness(physeq, "Health_status", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas+
  geom_boxplot(data=p_alpha_meas$data, aes(x=Health_status, y=value, fill=Health_status), fatten=NULL)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), 
               width= 0.75, linetype= "solid")+
  scale_fill_manual(values= c("#0077b6", "#ae2012"))+
  geom_jitter(width=0, size= 3.5, color= "black")+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  labs(title= "Chao1", y= "Diversity Index", x="Group")+
  scale_x_discrete(limits= c("Diseased", "Healthy"))

res.aov <- aov(value ~ Health_status, data= p_alpha_meas$data)  
summary(res.aov)

#Alpha diversity surface microbiomes HvD 
alpha_meas = c("ACE")
p_alpha_meas <- plot_richness(physeq, "Health_status", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas+
  geom_boxplot(data=p_alpha_meas$data, aes(x=Health_status, y=value, fill=Health_status), fatten=NULL)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), 
               width= 0.75, linetype= "solid")+
  scale_fill_manual(values= c("#0077b6", "#ae2012"))+
  geom_jitter(width=0, size= 3.5, color= "black")+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  labs(title= "Phylogenetic Diversity", y= "Diversity Index", x="Group")+
  scale_x_discrete(limits= c("Diseased", "Healthy"))

res.aov <- aov(value ~ Health_status, data= p_alpha_meas$data)  
summary(res.aov)

#NMDS on Weighted Unifrac distances
ordu_uni = ordinate(physeq, "NMDS", "unifrac", weighted=TRUE)
p_uni <- plot_ordination(physeq, ordu_uni, color="Health_status") +
  scale_color_manual(values= c("#0077b6", "#ae2012"))+
  geom_point(size=5) +
  ggtitle("MDS on weighted-UniFrac distance")
ggsave("16S_UNI_WEIGHTED_Silva_Status_C1.png", width = 8, height = 8, dpi = 600)

#add ellipses
p_uni + 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()
ggsave("16S_UNI_WEIGHTED_Silva_Status_ellipses_C1.png", width = 8, height = 8, dpi = 600)


##make barchart for most abundant taxa using 16S_OTU_TAX_merged_phylum
library(ggplot2)
library(dplyr)

data <- read.csv("Phylumabundance_surface.csv")


data_msd <- data %>% 
  group_by(Group, Taxa)%>%
  summarise_at(vars(Abundance),
               list( mean = mean,
                     sd=sd)) %>%
  as.data.frame()

#barplot per sample                                       
ggplot(data, aes(x= Sample, y= Abundance, fill=Taxa))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#ff7d00", "#ffc4d6", "#a4161a", "#31572c", "#74c69d", "#f1c453", "#a5a58d", 
                                        "#ff595e", "#4361ee", "#003566"))
                                        
## bargraph of species for surface microbiomes 
data <- read.csv("species_surface.csv")
data_msd <- data %>% 
  group_by(Group, Taxa)%>%
  summarise_at(vars(Abundance), list( mean = mean,
                                      sd=sd)) %>%
  as.data.frame()

#barplot per sample
ggplot(data, aes(x= Sample, y= Abundance, fill=Taxa))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#a4161a", "#613f75", "#c9184a", "#3f88c5", "#b5838d", "#ff7d00", "#ffc4d6", "#31572c", "#ff99c8",
                                        "#5c677d", "#74c69d", "#7f5539", "#cce6f4", "#a5a58d", "#ff595e", "#a98467", 
                                        "#6c584c","#de639a", "#40916c", "#fb8b24", "#168aad", "#e07a5f","#7400b8", "#e09f3e",
                                        "#3c6e71","#4361ee", "#9c89b8", "#acd8aa", "#003566"))
                                        

############################ C6 ######################################################

physeq <- physeq_C6

#get the data frame from the phyloseq object
pd <- psmelt(physeq)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)

# Changing order of levels
pd$Sample_type <- factor(pd$Sample_type,      # Reordering group factor levels
                         levels = c("lesion surface", "lesion body wall", "seawater swab"))
sample_data(physeq)$Sample_type <- factor(sample_data(physeq)$Sample_type, 
                                          levels = c("lesion surface", "lesion body wall", "seawater swab"), ordered = TRUE)
levels(sample_data(physeq)$Sample_type)


alpha_meas = c("Observed")
p_alpha_meas <- plot_richness(physeq, "Sample_type", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas+
  geom_boxplot(data=p_alpha_meas$data, aes(x=Sample_type, y=value, fill=Sample_type), fatten=NULL)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), 
               width= 0.75, linetype= "solid")+
  scale_fill_manual(values= c("#588157", "#ae2012", "#ff9f1c"))+
  geom_jitter(width=0, size= 3.5, color= "black")+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  labs(title= "Phylogenetic Diversity", y= "Diversity Index", x="Group")+
  scale_x_discrete(limits= c("lesion body wall", "lesion surface"))

alpha_meas = c("Chao1")
p_alpha_meas <- plot_richness(physeq, "Sample_type", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas+
  geom_boxplot(data=p_alpha_meas$data, aes(x=Sample_type, y=value, fill=Sample_type), fatten=NULL)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), 
               width= 0.75, linetype= "solid")+
  scale_fill_manual(values= c("#588157", "#ae2012", "#ff9f1c"))+
  geom_jitter(width=0, size= 3.5, color= "black")+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  labs(title= "Phylogenetic Diversity", y= "Diversity Index", x="Group")+
  scale_x_discrete(limits= c("lesion body wall", "lesion surface"))

alpha_meas = c("ACE")
p_alpha_meas <- plot_richness(physeq, "Sample_type", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas+
  geom_boxplot(data=p_alpha_meas$data, aes(x=Sample_type, y=value, fill=Sample_type), fatten=NULL)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), 
               width= 0.75, linetype= "solid")+
  scale_fill_manual(values= c("#588157", "#ae2012", "#ff9f1c"))+
  geom_jitter(width=0, size= 3.5, color= "black")+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  labs(title= "Phylogenetic Diversity", y= "Diversity Index", x="Group")+
  scale_x_discrete(limits= c("lesion body wall", "lesion surface"))

#NMDS on Weighted Unifrac distances
ordu_uni = ordinate(physeq, "NMDS", "unifrac", weighted=TRUE)
p_uni <- plot_ordination(physeq, ordu_uni, color="Health_status", shape = "Sample_type") +
  geom_point(size=5, alpha=0.75) +
  scale_colour_brewer(type="qual", palette="Set1")+ 
  ggtitle("MDS on weighted-UniFrac distance")
ggsave("16S_UNI_WEIGHTED_Silva_Status_C6.png", width = 8, height = 8, dpi = 600)

#add ellipses
p_uni + 
  stat_ellipse(type = "norm", linetype = 1) +
  #stat_ellipse(type = "t") +
  theme_bw()
ggsave("16S_UNI_WEIGHTED_Silva_Status_ellipses_C6.png", width = 8, height = 8, dpi = 600)

#NMDS on Weighted Unifrac distances
ordu_uni = ordinate(physeq, "NMDS", "unifrac", weighted=TRUE)
p_uni <- plot_ordination(physeq, ordu_uni, color="Sample_type") +
  scale_color_manual(values= c("#588157", "#ae2012"))+
  geom_point(size=5) +
  ggtitle("MDS on weighted-UniFrac distance")
ggsave("16S_UNI_WEIGHTED_Silva_Status_C1.png", width = 8, height = 8, dpi = 600)

#add ellipses
p_uni + 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()

##make barchart for most abundant taxa using 16S_OTU_TAX_merged_phylum
library(ggplot2)
library(dplyr)

data <- read.csv("phylaabundance_LBWLS.csv")


data_msd <- data %>% 
  group_by(Group, Taxa)%>%
  summarise_at(vars(Abundance),
               list( mean = mean,
                     sd=sd)) %>%
  as.data.frame()

#barplot per sample                                       
ggplot(data, aes(x= Sample, y= Abundance, fill=Taxa))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#ff7d00", "#ffc4d6", "#a4161a", "#31572c", "#74c69d", "#f1c453", "#a5a58d", 
                                        "#ff595e", "#4361ee", "#003566"))
                                        

## bargraph of species for surface microbiomes 
data <- read.csv("species_LBWvLS_nomito.csv")
data_msd <- data %>% 
  group_by(Group, Taxa)%>%
  summarise_at(vars(Abundance), list( mean = mean,
                                      sd=sd)) %>%
  as.data.frame()

#barplot per sample
ggplot(data, aes(x= Sample, y= Abundance, fill=Taxa))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#a4161a", "#613f75", "#c9184a", "#3f88c5", "#b5838d", "#ff7d00", "#ffc4d6", "#31572c", "#ff99c8",
                                        "#5c677d", "#74c69d", "#7f5539", "#cce6f4", "#a5a58d", "#ff595e", "#a98467", 
                                        "#de639a", "#40916c", "#fb8b24", "#168aad", "#e07a5f","#7400b8", "#e09f3e",
                                        "#3c6e71","#4361ee", "#9c89b8", "#acd8aa", "#003566"))
                                      



####################### C8 ####################################################

physeq <- physeq_C8

#get the data frame from the phyloseq object
pd <- psmelt(physeq)

#Count how many Phyla are there in your samples
HowManyPhyla <- length(unique(unlist(pd[,c("Phylum")])))

# Build a colour palette with number of colours as many as 
# the Phyla in your samples by interpolating the palette "Dark2".
getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
PhylaPalette = getPalette(HowManyPhyla)

# Changing order of levels
pd$Health_status <- factor(pd$Health_status,      # Reordering group factor levels
                           levels = c("Healthy", "Diseased"))
sample_data(physeq)$Health_status <- factor(sample_data(physeq)$Health_status, 
                                            levels = c("Healthy", "Diseased"), ordered = TRUE)
levels(sample_data(physeq)$Health_status)

alpha_meas = c("Observed")
p_alpha_meas <- plot_richness(physeq, "Group", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas+
  geom_boxplot(data=p_alpha_meas$data, aes(x=Group, y=value, fill=Group), fatten=NULL)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), 
               width= 0.75, linetype= "solid")+
  scale_fill_manual(values= c("#8f2d56", "#ffadad", "#89c2d9", "#014f86", "#ef233c"))+
  geom_jitter(width=0, size= 3.5, color= "black")+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  labs(title= "Phylogenetic Diversity", y= "Diversity Index", x="Group")+
  scale_x_discrete(limits= c("lesion body wall", "diseased body wall", "diseased coelomic fluid", "healthy body wall", "healthy coelomic fluid"))

res.aov <- aov(value ~ Group, data= p_alpha_meas$data)  
summary(res.aov)

alpha_meas = c("Chao1")
p_alpha_meas <- plot_richness(physeq, "Group", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas+
  geom_boxplot(data=p_alpha_meas$data, aes(x=Group, y=value, fill=Group), fatten=NULL)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), 
               width= 0.75, linetype= "solid")+
  scale_fill_manual(values= c("#8f2d56", "#ffadad", "#89c2d9", "#014f86", "#ef233c"))+
  geom_jitter(width=0, size= 3.5, color= "black")+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  labs(title= "Phylogenetic Diversity", y= "Diversity Index", x="Group")+
  scale_x_discrete(limits= c("lesion body wall", "diseased body wall", "diseased coelomic fluid", "healthy body wall", "healthy coelomic fluid"))

alpha_meas = c("ACE")
p_alpha_meas <- plot_richness(physeq, "Group", measures=alpha_meas)
#add a ggplot2 box plot layer to the previous plot
p_alpha_meas+
  geom_boxplot(data=p_alpha_meas$data, aes(x=Group, y=value, fill=Group), fatten=NULL)+
  stat_summary(fun = mean, geom = "errorbar", aes(ymax = after_stat(y), ymin = after_stat(y)), 
               width= 0.75, linetype= "solid")+
  scale_fill_manual(values= c("#8f2d56", "#ffadad", "#89c2d9", "#014f86", "#ef233c"))+
  geom_jitter(width=0, size= 3.5, color= "black")+
  theme(plot.title = element_text(hjust = 0.5, face= "bold"))+
  labs(title= "Phylogenetic Diversity", y= "Diversity Index", x="Group")+
  scale_x_discrete(limits= c("lesion body wall", "diseased body wall", "diseased coelomic fluid", "healthy body wall", "healthy coelomic fluid"))

#NMDS on Weighted Unifrac distances
ordu_uni = ordinate(physeq, "NMDS", "unifrac", weighted=TRUE)
p_uni <- plot_ordination(physeq, ordu_uni, color="Group")+ 
  scale_color_manual(values= c("#8f2d56", "#ffadad", "#89c2d9", "#014f86", "#ef233c"))+
  geom_point(size=4) +
  ggtitle("MDS on weighted-UniFrac distance")

#add ellipses
p_uni + 
  stat_ellipse(type = "norm", linetype = 2) +
  stat_ellipse(type = "t") +
  theme_bw()

#PERMANOVA in phyloseq
metadata_permanova <- as(sample_data(physeq), "data.frame")
permanova.Status <- adonis2(distance(physeq, method="bray") ~ Group, data = metadata_permanova)
permanova.Status


##make barchart for most abundant taxa using 16S_OTU_TAX_merged_phylum
library(ggplot2)
library(dplyr)

data <- read.csv("Phylumabundance_tissue.csv")


data_msd <- data %>% 
  group_by(Group, Taxa)%>%
  summarise_at(vars(Abundance),
               list( mean = mean,
                     sd=sd)) %>%
  as.data.frame()

#barplot per sample                                       
ggplot(data, aes(x= Sample, y= Abundance, fill=Taxa))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#ff7d00", "#ffc4d6", "#a4161a", "#31572c", "#74c69d", "#f1c453", "#a5a58d", 
                                        "#ff595e", "#4361ee", "#003566"))+
                                          scale_x_discrete(limits= c("D1LBW", "D2aLBW", "D2bLBW", "D3LBW", "D1BW", 
                                                                     "D2BW", "D3BW", "D1CF", "D2CF", "D3CF", "H1BW", 
                                                                     "H2BW", "H3BW", "H1CF", "H2CF", "H3CF"))


## bargraph of species for surface microbiomes 
data <- read.csv("species_tissue_nomito.csv")
data_msd <- data %>% 
  group_by(Group, Taxa)%>%
  summarise_at(vars(Abundance), list( mean = mean,
                                      sd=sd)) %>%
  as.data.frame()

#barplot per sample
ggplot(data, aes(x= Sample, y= Abundance, fill=Taxa))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))+
  scale_fill_manual(values = c("#a4161a", "#613f75", "#c9184a", "#3f88c5", "#b5838d", "#ff7d00", "#ffc4d6", "#31572c", "#ff99c8",
                                        "#5c677d", "#74c69d", "#7f5539", "#cce6f4", "#a5a58d", "#ff595e", "#a98467", 
                                        "#de639a", "#40916c", "#fb8b24", "#168aad", "#e07a5f","#7400b8", "#e09f3e",
                                        "#3c6e71","#4361ee", "#9c89b8", "#acd8aa", "#003566"))+
                                          scale_x_discrete(limits= c("D1LBW", "D2aLBW", "D2bLBW", "D3LBW", "D1BW", 
                                                                     "D2BW", "D3BW", "D1CF", "D2CF", "D3CF", "H1BW", 
                                                                     "H2BW", "H3BW", "H1CF", "H2CF", "H3CF"))



