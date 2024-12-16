library(ALDEx2)
library(phyloseq)
library(ape) 
library(tidyverse)
library(picante)
library(ggplot2)
library(microbiome)

###############################################################################
# Step 1 - Load Phyloseq Object #########################################
###############################################################################

load("MS_phyloseq.RData")

###############################################################################
# Step 2 - Data Preparation (Binning) #########################################
###############################################################################

# 1. Bin phyloseq data based on age group #####################################

MSphylo_filt@sam_data$agegroup <- cut(MSphylo_filt@sam_data$age, 
                                      breaks = c(25, 55, Inf), 
                                      labels = c("young", "old"),
                                      right = FALSE)
MSphylo_filt@sam_data$agegroup

sample_data(MSphylo_filt)$age_disease <- ifelse(
  sample_data(MSphylo_filt)$agegroup == "old" & sample_data(MSphylo_filt)$disease_course %in% c('PPMS', 'SPMS'), "old, PMS",
  ifelse(
    sample_data(MSphylo_filt)$agegroup == "young" & sample_data(MSphylo_filt)$disease_course %in% c('PPMS', 'SPMS'), "young, PMS",
    ifelse(
      sample_data(MSphylo_filt)$agegroup == "old" & sample_data(MSphylo_filt)$disease_course == 'Control', "old, healthy",
      ifelse(
        sample_data(MSphylo_filt)$agegroup == "young" & sample_data(MSphylo_filt)$disease_course == 'Control', "young, healthy",
        NA ))))
get_variable(MSphylo_filt, c("age_disease"))

# 2. Make a column for age, disease, and smoking status #######################

sample_data(MSphylo_filt)$age_disease_smoking <- ifelse(
  sample_data(MSphylo_filt)$age_disease == "old, healthy" & sample_data(MSphylo_filt)$smoke %in% c('formersmoker', 'smoker'), "old, healthy, smoker", 
  ifelse(
    sample_data(MSphylo_filt)$age_disease == "young, healthy" & sample_data(MSphylo_filt)$smoke %in% c('formersmoker', 'smoker'), "young, healthy, smoker",
    ifelse(
      sample_data(MSphylo_filt)$age_disease == "old, healthy" & sample_data(MSphylo_filt)$smoke == 'nonsmoker', "old, healthy, nonsmoker",
      ifelse(
        sample_data(MSphylo_filt)$age_disease == "young, healthy" & sample_data(MSphylo_filt)$smoke == 'nonsmoker', "young, healthy, nonsmoker",
        ifelse(
          sample_data(MSphylo_filt)$age_disease == "old, PMS" & sample_data(MSphylo_filt)$smoke == 'nonsmoker', "old, PMS, nonsmoker",
          ifelse(
            sample_data(MSphylo_filt)$age_disease == "young, PMS" & sample_data(MSphylo_filt)$smoke == 'nonsmoker', "young, PMS, nonsmoker",
            ifelse(
              sample_data(MSphylo_filt)$age_disease == "old, PMS" & sample_data(MSphylo_filt)$smoke %in% c('formersmoker', 'smoker'), "old, PMS, smoker",
              ifelse(
                sample_data(MSphylo_filt)$age_disease == "young, PMS" & sample_data(MSphylo_filt)$smoke %in% c('formersmoker', 'smoker'), "young, PMS, smoker", NA ))))))))
get_variable(MSphylo_filt, c("age_disease_smoking"))

#metadata = sample_data(MSphylo_filt)
#metadata

################################################################################
# Step 2 - Filter Data #########################################################
################################################################################

# Young Nonsmokers Healthy v. pMS
subset_young_nonsmoker <- subset_samples(MSphylo_filt, age_disease_smoking == "young, healthy, nonsmoker" | 
                                   age_disease_smoking == "young, PMS, nonsmoker")
subset_young_nonsmoker

# Young Smokers Healthy v. pMS
subset_young_smoker <- subset_samples(MSphylo_filt, age_disease_smoking == "young, healthy, smoker" | 
                                           age_disease_smoking == "young, PMS, smoker")
subset_young_smoker

# Older Nonsmokers Healthy v. pMS
subset_old_nonsmoker <- subset_samples(MSphylo_filt, age_disease_smoking == "old, healthy, nonsmoker" | 
                                         age_disease_smoking == "old, PMS, nonsmoker")
subset_old_nonsmoker

# Older Smokers Healthy v. pMS
subset_old_smoker <- subset_samples(MSphylo_filt, age_disease_smoking == "old, healthy, smoker" | 
                                         age_disease_smoking == "old, PMS, smoker")
subset_old_smoker

# aggregation at the family level
family_young_nonsmoker <- tax_glom(subset_young_nonsmoker,'Family')
ntaxa(subset_young_nonsmoker); ntaxa(family_young_nonsmoker)

family_young_smoker <- tax_glom(subset_young_smoker,'Family')
ntaxa(subset_young_smoker); ntaxa(family_young_smoker)

family_old_smoker <- tax_glom(subset_old_smoker,'Family')
ntaxa(subset_old_smoker); ntaxa(family_old_smoker)

family_old_nonsmoker <- tax_glom(subset_old_nonsmoker,'Family')
ntaxa(subset_old_nonsmoker); ntaxa(family_old_nonsmoker)

# aggregation at the genus level
genus_young_nonsmoker <- tax_glom(subset_young_nonsmoker,'Genus')
ntaxa(subset_young_nonsmoker); ntaxa(genus_young_nonsmoker)

genus_young_smoker <- tax_glom(subset_young_smoker,'Genus')
ntaxa(subset_young_smoker); ntaxa(genus_young_nonsmoker)

genus_old_smoker <- tax_glom(subset_old_smoker,'Genus')
ntaxa(subset_old_smoker); ntaxa(genus_young_smoker)

genus_old_nonsmoker <- tax_glom(subset_old_nonsmoker,'Genus')
ntaxa(subset_old_nonsmoker); ntaxa(genus_young_nonsmoker)

# filtering out ASVs (only including ASVs that are more than 0.01%)
calculate_relative_abundance <- function(x) x / sum(x)

######## FAMILY LEVEL ##########################################################
family_young_smoker_counts <- taxa_sums(family_young_smoker) # Total sum for that taxa
family_young_smoker_counts_abundance <- calculate_relative_abundance(family_young_smoker_counts) # overall proportion of each bug
family_young_smoker_abundance <- family_young_smoker_counts_abundance > 0.001 # is each bug above the threshold? TRUE if so.
family_young_smoker <- prune_taxa(family_young_smoker_abundance, family_young_smoker) # Take only bugs above threshold
family_young_smoker
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 29 taxa and 5 samples ]
#sample_data() Sample Data:       [ 5 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 29 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 29 tips and 28 internal nodes ]
family_young_nonsmoker_counts <- taxa_sums(family_young_nonsmoker) # Total sum for that taxa
family_young_nonsmoker_counts_abundance <- calculate_relative_abundance(family_young_nonsmoker_counts) # overall proportion of each bug
family_young_nonsmoker_abundance <- family_young_nonsmoker_counts_abundance > 0.001 # is each bug above the threshold? TRUE if so.
family_young_nonsmoker <- prune_taxa(family_young_nonsmoker_abundance, family_young_nonsmoker) # Take only bugs above threshold
family_young_nonsmoker
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 32 taxa and 3 samples ]
#sample_data() Sample Data:       [ 3 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 32 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 32 tips and 31 internal nodes ]
family_old_smoker_counts <- taxa_sums(family_old_smoker)
family_old_smoker_counts_abundance <- calculate_relative_abundance(family_old_smoker_counts) 
family_old_smoker_abundance <- family_old_smoker_counts_abundance > 0.001
family_old_smoker <- prune_taxa(family_old_smoker_abundance, family_old_smoker)
family_old_smoker
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 37 taxa and 31 samples ]
#sample_data() Sample Data:       [ 31 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 37 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 37 tips and 36 internal nodes ]
family_old_nonsmoker_counts <- taxa_sums(family_old_nonsmoker)
family_old_nonsmoker_counts_abundance <- calculate_relative_abundance(family_old_nonsmoker_counts) 
family_old_nonsmoker_abundance <- family_old_nonsmoker_counts_abundance > 0.001
family_old_nonsmoker <- prune_taxa(family_old_nonsmoker_abundance, family_old_nonsmoker)
family_old_nonsmoker
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 40 taxa and 71 samples ]
#sample_data() Sample Data:       [ 71 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 40 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 40 tips and 39 internal nodes ]
########## GENUS LEVEL #########################################################
genus_young_smoker_counts <- taxa_sums(genus_young_smoker) # Total sum for that taxa
genus_young_smoker_counts_abundance <- calculate_relative_abundance(genus_young_smoker_counts) # overall proportion of each bug
genus_young_smoker_abundance <- genus_young_smoker_counts_abundance > 0.001 # is each bug above the threshold? TRUE if so.
genus_young_smoker <- prune_taxa(genus_young_smoker_abundance, genus_young_smoker) # Take only bugs above threshold
genus_young_smoker
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 66 taxa and 5 samples ]
#sample_data() Sample Data:       [ 5 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 66 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 66 tips and 65 internal nodes ]
genus_young_nonsmoker_counts <- taxa_sums(genus_young_nonsmoker) # Total sum for that taxa
genus_young_nonsmoker_counts_abundance <- calculate_relative_abundance(genus_young_nonsmoker_counts) # overall proportion of each bug
genus_young_nonsmoker_abundance <- genus_young_nonsmoker_counts_abundance > 0.001 # is each bug above the threshold? TRUE if so.
genus_young_nonsmoker <- prune_taxa(genus_young_nonsmoker_abundance, genus_young_nonsmoker) # Take only bugs above threshold
genus_young_nonsmoker
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 68 taxa and 3 samples ]
#sample_data() Sample Data:       [ 3 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 68 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 68 tips and 67 internal nodes ]
genus_old_smoker_counts <- taxa_sums(genus_old_smoker)
genus_old_smoker_counts_abundance <- calculate_relative_abundance(genus_old_smoker_counts) 
genus_old_smoker_abundance <- genus_old_smoker_counts_abundance > 0.001
genus_old_smoker <- prune_taxa(genus_old_smoker_abundance, genus_old_smoker)
genus_old_smoker
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 85 taxa and 31 samples ]
#sample_data() Sample Data:       [ 31 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 85 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 85 tips and 84 internal nodes ]
genus_old_nonsmoker_counts <- taxa_sums(genus_old_nonsmoker)
genus_old_nonsmoker_counts_abundance <- calculate_relative_abundance(genus_old_nonsmoker_counts) 
genus_old_nonsmoker_abundance <- genus_old_nonsmoker_counts_abundance > 0.001
genus_old_nonsmoker <- prune_taxa(genus_old_nonsmoker_abundance, genus_old_nonsmoker)
genus_old_nonsmoker
#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 92 taxa and 71 samples ]
#sample_data() Sample Data:       [ 71 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 92 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 92 tips and 91 internal nodes ]



###############################################################################
# Step 4 - ALDeX2 #############################################################
###############################################################################

# ALDeX2: testing differential abundance across multiple conditions or groups and 
######### adjusts for within-sample variation.

#set.seed(124) #for reproducibility

#################################################################################
## FAMILY LEVEL #################################################################
#################################################################################

# Old Smokers Healthy v. pMS ################################################
s = family_old_smoker@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = family_old_smoker@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df= aldex.glm(x) #old_pms_family
colnames(df)
# "age_disease_smokingold, PMS, smoker:Est" = log_fold_change
# "age_disease_smokingold, PMS, smoker:pval" = p-value
# negative means that the ASVs are down-regulated

log_fold_change <- df$`age_disease_smokingold, PMS, smoker:Est`
pvals <- df$`age_disease_smokingold, PMS, smoker:pval`
length(log_fold_change)
length(pvals)   

valid_idx <- !is.na(log_fold_change) & !is.na(pvals)

log_fold_change <- log_fold_change[valid_idx]
pvals <- pvals[valid_idx]

volcano_data <- data.frame(
  logFC = log_fold_change,                 # Log fold change
  pval = pvals,                            # p-values
  logP = -log10(pvals),                    # -log10 of p-values (for better visualization)
  Significant = ifelse(pvals < 0.05 & abs(log_fold_change) > 1, "Yes", "No") # Significance threshold
)

ggplot(volcano_data, aes(x = logFC, y = logP, color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Family Level - Older Smokers Healthy v. pMS", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# SIGNIFICANCE = several families are decreased and increased

# Old Nonsmokers Healthy v. pMS ################################################
s = family_old_nonsmoker@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = family_old_nonsmoker@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df = aldex.glm(x) # young_pms_family
colnames(df)

log_fold_change <- df$`age_disease_smokingold, PMS, nonsmoker:Est`
pvals <- df$`age_disease_smokingold, PMS, nonsmoker:pval`

valid_idx <- !is.na(log_fold_change) & !is.na(pvals)

log_fold_change <- log_fold_change[valid_idx]
pvals <- pvals[valid_idx]

volcano_data <- data.frame(
  logFC = log_fold_change,                 # Log fold change
  pval = pvals,                            # p-values
  logP = -log10(pvals),                    # -log10 of p-values (for better visualization)
  Significant = ifelse(pvals < 0.05 & abs(log_fold_change) > 1, "Yes", "No") # Significance threshold
)

ggplot(volcano_data, aes(x = logFC, y = logP, color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Family Level - Older Nonsmokers Healthy v. pMS", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# SIGNIFICANCE = 1 family decrease

###### Young Nonsmokers Healthy v. pMS #############################################
s = family_young_nonsmoker@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = family_young_nonsmoker@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df = aldex.glm(x) 
colnames(df)

log_fold_change <- df$`age_disease_smokingyoung, PMS, nonsmoker:Est`
pvals <- df$`age_disease_smokingyoung, PMS, nonsmoker:pval`

valid_idx <- !is.na(log_fold_change) & !is.na(pvals)

log_fold_change <- log_fold_change[valid_idx]
pvals <- pvals[valid_idx]

volcano_data <- data.frame(
  logFC = log_fold_change,                 # Log fold change
  pval = pvals,                            # p-values
  logP = -log10(pvals),                    # -log10 of p-values (for better visualization)
  Significant = ifelse(pvals < 0.05 & abs(log_fold_change) > 1, "Yes", "No") # Significance threshold
)

ggplot(volcano_data, aes(x = logFC, y = logP, color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Family Level - Younger Nonsmokers: Healthy v. pMS ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# NOTE. there seems to be TWO familiies of relevance

# Young Smokers Healthy v. pMS ############
s = family_young_smoker@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = family_young_smoker@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df = aldex.glm(x) 
colnames(df)

log_fold_change <- df$`age_disease_smokingyoung, PMS, smoker:Est`
pvals <- df$`age_disease_smokingyoung, PMS, smoker:pval`

valid_idx <- !is.na(log_fold_change) & !is.na(pvals)

log_fold_change <- log_fold_change[valid_idx]
pvals <- pvals[valid_idx]

volcano_data <- data.frame(
  logFC = log_fold_change,                 # Log fold change
  pval = pvals,                            # p-values
  logP = -log10(pvals),                    # -log10 of p-values (for better visualization)
  Significant = ifelse(pvals < 0.05 & abs(log_fold_change) > 1, "Yes", "No") # Significance threshold
)

ggplot(volcano_data, aes(x = logFC, y = logP, color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Family Level - Younger Smokers: Healthy v. pMS ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")
# NOTE. there seems to be 3 family of relevance

################################################################################
## GENUS LEVEL #################################################################
################################################################################

###############################################################################
################## Older Smoker Healthy v. pMS #############################
###############################################################################
s = genus_old_smoker@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = genus_old_smoker@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df = aldex.glm(x) 
colnames(df)

log_fold_change <- df$`age_disease_smokingold, PMS, smoker:Est`
pvals <- df$`age_disease_smokingold, PMS, smoker:pval`

valid_idx <- !is.na(log_fold_change) & !is.na(pvals)

log_fold_change <- log_fold_change[valid_idx]
pvals <- pvals[valid_idx]

volcano_data <- data.frame(
  logFC = log_fold_change,                 # Log fold change
  pval = pvals,                            # p-values
  logP = -log10(pvals),                    # -log10 of p-values (for better visualization)
  Significant = ifelse(pvals < 0.05 & abs(log_fold_change) > 1, "Yes", "No") # Significance threshold
)

ggplot(volcano_data, aes(x = logFC, y = logP, color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Genus Level - Older Smokers: Healthy v. pMS ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

###############################################################################
################## Older Nonsmoker Healthy v. pMS #############################
###############################################################################
s = genus_old_nonsmoker@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = genus_old_nonsmoker@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df = aldex.glm(x) # young_pms_genus
colnames(df)

log_fold_change <- df$`age_disease_smokingold, PMS, nonsmoker:Est`
pvals <- df$`age_disease_smokingold, PMS, nonsmoker:pval`

valid_idx <- !is.na(log_fold_change) & !is.na(pvals)

log_fold_change <- log_fold_change[valid_idx]
pvals <- pvals[valid_idx]

volcano_data <- data.frame(
  logFC = log_fold_change,                 # Log fold change
  pval = pvals,                            # p-values
  logP = -log10(pvals),                    # -log10 of p-values (for better visualization)
  Significant = ifelse(pvals < 0.05 & abs(log_fold_change) > 1, "Yes", "No") # Significance threshold
)

ggplot(volcano_data, aes(x = logFC, y = logP, color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Genus Level - Older Nonsmokers: Healthy v. pMS", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

###############################################################################
################## Young Smoker Healthy v. pMS ################################
###############################################################################
s = genus_young_smoker@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = genus_young_smoker@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df = aldex.glm(x) # young_healthy_family
colnames(df)

log_fold_change <- df$`age_disease_smokingyoung, PMS, smoker:Est`
pvals <- df$`age_disease_smokingyoung, PMS, smoker:pval`

valid_idx <- !is.na(log_fold_change) & !is.na(pvals)

log_fold_change <- log_fold_change[valid_idx]
pvals <- pvals[valid_idx]

volcano_data <- data.frame(
  logFC = log_fold_change,                 # Log fold change
  pval = pvals,                            # p-values
  logP = -log10(pvals),                    # -log10 of p-values (for better visualization)
  Significant = ifelse(pvals < 0.05 & abs(log_fold_change) > 1, "Yes", "No") # Significance threshold
)

ggplot(volcano_data, aes(x = logFC, y = logP, color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Genus Level - Younger Smokers: Healthy v. pMS ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

###############################################################################
################## Young Nonsmoker Healthy v. pMS #############################
###############################################################################
s = genus_young_nonsmoker@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = genus_young_nonsmoker@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df = aldex.glm(x) # young_healthy_family
colnames(df)

log_fold_change <- df$`age_disease_smokingyoung, PMS, nonsmoker:Est`
pvals <- df$`age_disease_smokingyoung, PMS, nonsmoker:pval`

valid_idx <- !is.na(log_fold_change) & !is.na(pvals)

log_fold_change <- log_fold_change[valid_idx]
pvals <- pvals[valid_idx]

volcano_data <- data.frame(
  logFC = log_fold_change,                 # Log fold change
  pval = pvals,                            # p-values
  logP = -log10(pvals),                    # -log10 of p-values (for better visualization)
  Significant = ifelse(pvals < 0.05 & abs(log_fold_change) > 1, "Yes", "No") # Significance threshold
)

ggplot(volcano_data, aes(x = logFC, y = logP, color = Significant)) +
  geom_point(alpha = 0.6, size = 2) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(title = "Genus Level - Younger Nonsmokers: Healthy v. pMS ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")
