# MICB 475 - TEAM 10 - Burak's Processing & Diversity Metrics File

###### Load in necessary libraries ############################################
library(phyloseq)
library(ape) 
library(tidyverse)
library(picante)
library(ggplot2)
library(vegan)     # for distance and PCoA
library(ggplot2)   # for ggplot2
library(ggalt)     # for geom_encircle

###### Load in RData ###########################################################
meta_new <- "ms_metadata.csv"
meta <- read.csv(meta_new)
load("MS_rare.RData") # the data file is called "MS_rare"
ls()
otu_table(MS_rare)
sample_data(MS_rare)
tax_table(MS_rare)
phy_tree(MS_rare)

###### Accessing Rarefied Phyloseq Object ######################################
sample_variables(MS_rare)
get_variable(MS_rare, c("age", "disease_course")) # note that disease course has SPMS or PPMS
sample_data_df <- as.data.frame(sample_data(MS_rare))

####### Processing the Data ####################################################
MSrare_filt <- subset_taxa(MS_rare, Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria") # getting rid of all the species that are not bacterial
MSrare_filt_nolow <- filter_taxa(MSrare_filt, function(x) sum(x)>5, prune = TRUE) # remove ASVs with counts lower than 5
MSrare_filt_nolow_samps <- prune_samples(sample_sums(MSrare_filt_nolow)>100, MSrare_filt_nolow) # remove samples with reads below 100
MSrare_almost <- subset_samples(MSrare_filt_nolow_samps, !is.na(disease_course)) # no disease_course NAs
MSrare_final <- subset_samples(MSrare_almost, !is.na(age)) # no age NAs
# no rarefaction needed as this has already been conducted

####### Grouping the Data ######################################################
# ifelse(condition, value_if_true, value_if_false)
# basically I want to create a new column where the samples are given a new identity
# young_healthy, young_pMS, old_healthy, old_pMS
# young is less than 55 years of age
# older is 55 and above
# young is Control
# pMS is SPMS or PPMS
sample_data(MSrare_final)$age_disease <- ifelse(
  sample_data(MSrare_final)$disease_course %in% c("SPMS", "PPMS"),
  ifelse(sample_data(MSrare_final)$age < 55, "young_pMS", "older_pMS"),
  ifelse(sample_data(MSrare_final)$age < 55, "young_healthy", "older_healthy")
)
####### Diversity Metrics ######################################################
# try all the beta diversity metrics
# Bray-Curtis - takes abundance (evenness) into account
bc_dm <- distance(MSrare_final, method="bray")
pcoa_bc <- ordinate(MSrare_final, method="PCoA", distance=bc_dm)
plot_ordination(MSrare_final, pcoa_bc, color = "age_disease")
gg_pcoa <- plot_ordination(MSrare_final, pcoa_bc, color = "age_disease") +
  labs(pch="Age Disease", col = "Groups") +
  stat_ellipse(type = "norm") +
  labs(title="Beta Diversity of pMS v. Healthy Groups (Bray-Curtis Distance)") +
  theme_classic() +
  scale_color_manual(values = c("older_healthy" = "blue", 
                                "young_healthy"= "purple",
                                "older_pMS" = "gray",
                                "young_pMS" = "black"), 
                     labels = c("young_healthy" = "Younger Healthy",
                                "older_healthy" = "Older Healthy", 
                                "young_pMS" = "Younger pMS",
                                "older_pMS" = "Older pMS"))
gg_pcoa

# Jaccard - only cares about presence/absence
jc_dm <- distance(MSrare_final, method="jaccard")
pcoa_jc <- ordinate(MSrare_final, method="PCoA", distance=jc_dm)

p <- plot_ordination(MSrare_final, pcoa_jc, color = "age_disease")
p + 
   geom_encircle(aes(x = PC1, y = PC2, color = age_disease), 
                 data = p$data, 
                 expand = 0.1, 
                 show.legend = FALSE) +
   theme_minimal() +
   theme(legend.position = "top") +
   labs(title = "PCoA with Circles Around Groups")

# Unweighted UniFrac Distance - cares about presence/absence, and also relatedness (phylogenic distance)
uu_dm <- distance(MSrare_final, method="unifrac", weighted=FALSE)
pcoa_uu <- ordinate(MSrare_final, method="PCoA", distance=uu_dm)
plot_ordination(MSrare_final, pcoa_uu, color = "age_disease")

# Weighted UniFrac Distance - cares about abundance AND relatedness (phylogenic distance)
wu_dm <- distance(MSrare_final, method="unifrac", weighted=TRUE)
pcoa_wu <- ordinate(MSrare_final, method="PCoA", distance=wu_dm)
plot_ordination(MSrare_final, pcoa_wu, color = "age_disease")

# 1. Statistical Analyses
# 2. Write up the results




# Improve the plots with ggplot2

gg_pcoa <- plot_ordination(MS_rare, pcoa_bc, color = "body.site", shape="subject") +
  labs(pch="Subject #", col = "Body Site")
gg_pcoa

ggsave("plot_pcoa.png", gg_pcoa, height=4, width=5, dpi=300)

###############################################################################
# Identifying differentially present species across sample groups via LDEX2





