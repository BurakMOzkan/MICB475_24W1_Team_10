# analysis of microbiome composition  - differential abundance analysis
#package installation
install.packages("BiocManager")
BiocManager::install("ALDEx2")

#load ALDEX2 package
library(ALDEx2)
library(phyloseq)
library(ape) 
library(tidyverse)
library(picante)
library(ggplot2)
library(microbiome)

#load phyloseq object
load("C:/Users/burak/Desktop/Team10MS/MS_phyloseq.RData")

###############################################################################
# Step 1 - Data Preparation ###################################################

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

###############################################################################
# Step 2 - Filter Data ###################################################
# only include samples that are in the filtered metadata

#comparing old, PMS, smokers to nonsmokers
subset_old.pms <- subset_samples(MSphylo_filt, age_disease_smoking == "old, PMS, smoker" | 
                                   age_disease_smoking == "old, PMS, nonsmoker")
subset_old.pms

#comparing young, PMS, smokers to nonsmokers
subset_young.pms <- subset_samples(MSphylo_filt, age_disease_smoking == "young, PMS, smoker" | 
                                     age_disease_smoking == "young, PMS, nonsmoker")

#comparing old, healthy, smokers to nonsmokers
subset_old.healthy <- subset_samples(MSphylo_filt, age_disease_smoking == "old, healthy, smoker" | 
                                       age_disease_smoking == "old, healthy, nonsmoker")

#comparing young, healthy, smokers to nonsmokers
subset_young.healthy <- subset_samples(MSphylo_filt, age_disease_smoking == "young, healthy, smoker" | 
                                         age_disease_smoking == "young, healthy, nonsmoker")

#aggregate data at the family level to decrease amount of samples to analyse

# FAMILY LEVEL

family.old.pms <- tax_glom(subset_old.pms,'Family')
ntaxa(subset_old.pms); ntaxa(family.old.pms)

family.young.pms <- tax_glom(subset_young.pms,'Family')
ntaxa(subset_young.pms); ntaxa(family.young.pms)

family.old.healthy <- tax_glom(subset_old.healthy,'Family')
ntaxa(subset_old.healthy); ntaxa(family.old.healthy)

family.young.healthy <- tax_glom(subset_young.healthy,'Family')
ntaxa(subset_young.healthy); ntaxa(family.young.healthy)

#aggregate data at the genus level to decrease amount of samples to analyse

# GENUS LEVEL

genus.old.pms <- tax_glom(subset_old.pms,'Genus')
ntaxa(subset_old.pms); ntaxa(genus.old.pms)

genus.young.pms <- tax_glom(subset_young.pms,'Genus')
ntaxa(subset_young.pms); ntaxa(genus.young.pms)

genus.old.healthy <- tax_glom(subset_old.healthy,'Genus')
ntaxa(subset_old.healthy); ntaxa(genus.old.healthy)

genus.young.healthy <- tax_glom(subset_young.healthy,'Genus')
ntaxa(subset_young.healthy); ntaxa(genus.young.healthy)

#filtering out ASVs (only including ASVs that are more than 0.01%)

install.packages("microbiome")
library(microbiome)

# only include things that are at least 0.1% abundant (0.001) across all samples
calculate_relative_abundance <- function(x) x / sum(x)

#at the family level
family.old.pms.counts <- taxa_sums(family.old.pms) # Total sum for that taxa
family.old.pms.counts.abundance <- calculate_relative_abundance(family.old.pms.counts) # overall proportion of each bug
family.old.pms.abundant <- family.old.pms.counts.abundance > 0.001 # is each bug above the threshold? TRUE if so.
family.old.pms <- prune_taxa(family.old.pms.abundant, family.old.pms) # Take only bugs above threshold
family.old.pms 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 93 taxa and 34 samples ]
#sample_data() Sample Data:       [ 34 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 93 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 93 tips and 92 internal nodes ]

family.young.pms.counts <- taxa_sums(family.young.pms)
family.young.pms.counts.abundance <- calculate_relative_abundance(family.young.pms.counts) 
family.young.pms.abundant <- family.young.pms.counts.abundance > 0.001
family.young.pms <- prune_taxa(family.young.pms.abundant, family.young.pms)
family.young.pms 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 93 taxa and 8 samples ]
#sample_data() Sample Data:       [ 8 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 93 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 93 tips and 92 internal nodes ]

family.old.healthy.counts <- taxa_sums(family.old.healthy)
family.old.healthy.counts.abundance <- calculate_relative_abundance(family.old.healthy.counts) 
family.old.healthy.abundant <- family.old.healthy.counts.abundance > 0.001
family.old.healthy <- prune_taxa(family.old.healthy.abundant, family.old.healthy)
family.old.healthy 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 36 taxa and 68 samples ]
#sample_data() Sample Data:       [ 68 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 36 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 36 tips and 35 internal nodes ]

family.young.healthy.counts <- taxa_sums(family.young.healthy)
family.young.healthy.counts.abundance <- calculate_relative_abundance(family.young.healthy.counts) 
family.young.healthy.abundant <- family.young.healthy.counts.abundance > 0.001
family.young.healthy <- prune_taxa(family.young.healthy.abundant, family.young.healthy)
family.young.healthy

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 36 taxa and 35 samples ]
#sample_data() Sample Data:       [ 35 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 36 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 36 tips and 35 internal nodes ]

#at the genus level
genus.old.pms.counts <- taxa_sums(genus.old.pms) 
genus.old.pms.counts.abundance <- calculate_relative_abundance(genus.old.pms.counts) 
genus.old.pms.abundant <- genus.old.pms.counts.abundance > 0.001 
genus.old.pms <- prune_taxa(genus.old.pms.abundant, genus.old.pms) 
genus.old.pms 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 91 taxa and 34 samples ]
#sample_data() Sample Data:       [ 34 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 91 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 91 tips and 90 internal nodes ]

genus.young.pms.counts <- taxa_sums(genus.young.pms)
genus.young.pms.counts.abundance <- calculate_relative_abundance(genus.young.pms.counts) 
genus.young.pms.abundant <- genus.young.pms.counts.abundance > 0.001
genus.young.pms <- prune_taxa(genus.young.pms.abundant, genus.young.pms)
genus.young.pms 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 68 taxa and 8 samples ]
#sample_data() Sample Data:       [ 8 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 68 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 68 tips and 67 internal nodes ]

genus.old.healthy.counts <- taxa_sums(genus.old.healthy)
genus.old.healthy.counts.abundance <- calculate_relative_abundance(genus.old.healthy.counts) 
genus.old.healthy.abundant <- genus.old.healthy.counts.abundance > 0.001
genus.old.healthy <- prune_taxa(genus.old.healthy.abundant, genus.old.healthy)
genus.old.healthy 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 84 taxa and 68 samples ]
#sample_data() Sample Data:       [ 68 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 84 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 84 tips and 83 internal nodes ]

genus.young.healthy.counts <- taxa_sums(genus.young.healthy)
genus.young.healthy.counts.abundance <- calculate_relative_abundance(genus.young.healthy.counts) 
genus.young.healthy.abundant <- genus.young.healthy.counts.abundance > 0.001
genus.young.healthy <- prune_taxa(genus.young.healthy.abundant, genus.young.healthy)
genus.young.healthy

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 81 taxa and 35 samples ]
#sample_data() Sample Data:       [ 35 samples by 62 sample variables ]
#tax_table()   Taxonomy Table:    [ 81 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 81 tips and 80 internal nodes ]

###############################################################################
# Step 3 - ANCOM ###################################################

# ANCOM: non-parametric approach. relative abundance of each feature is 
######## constrained (i.e., the total of all taxa adds up to 1). comparing the 
######## relative abundance of each feature across different groups (conditions).
######## tests the null hypothesis that the abundance of a feature does not differ 
######## between groups. It uses a set of pairwise comparisons between features 
######## (in terms of their log-ratios) to detect whether there are significant differences in relative abundance.
######## a list of features with p-values indicating whether their relative abundances 
######## differ significantly across the groups being tested. Features with 
######## p-values below a certain threshold (typically 0.05) are considered 
######## differentially abundant.

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")
library(ANCOMBC)

###############################################################################
# FAMILY LEVEL ################################################################
###############################################################################

# old_pms smoking v. old_pms non-smoking

ancom.family.old.pms = ancombc(phyloseq = family.old.pms, # aggregated phyloseq
                               formula = 'age_disease_smoking', # The explanatory variable
                               p_adj_method = "fdr", # false discovery rate
                               prv_cut=0.10, # Max proportion of zeros allowed per taxon, get rid of high values of zeros
                               lib_cut = 1000, # Can filter out samples below minimum sequence depth here
                               group = 'age_disease_smoking', # If you're including structural zeros below, you need this
                               struc_zero = T) # If true, any taxa present in only 1 category of 'group' will automatically be significant
str(ancom.family.old.pms)
View(ancom.family.old.pms)

# young_pms smoking v. young_pms non-smoking ##################################

ancom.family.young.pms = ancombc(phyloseq = family.young.pms, 
                                 formula = 'age_disease_smoking', 
                                 p_adj_method = "fdr",
                                 prv_cut=0.10, 
                                 lib_cut = 1000, 
                                 group = 'age_disease_smoking', 
                                 struc_zero = T) 
str(ancom.family.young.pms)
View(ancom.family.young.pms)

# young_healthy smoking v. young_healthy non-smoking

ancom.family.young.healthy = ancombc(phyloseq = family.young.healthy,
                                     formula = 'age_disease_smoking',
                                     p_adj_method = "fdr",
                                     prv_cut=0.10, 
                                     lib_cut = 1000, 
                                     group = 'age_disease_smoking', 
                                     struc_zero = T) 

# old_healthy smoking v. old_healthy non-smoking

ancom.family.old.healthy = ancombc(phyloseq = family.old.healthy, 
                                   formula = 'age_disease_smoking', 
                                   p_adj_method = "fdr",
                                   prv_cut=0.10, 
                                   lib_cut = 1000,
                                   group = 'age_disease_smoking', 
                                   struc_zero = T)
str(ancom.family.old.healthy)
View(ancom.family.old.healthy)

###############################################################################
# GENUS LEVEL #################################################################
###############################################################################

# old_pms smoking v. old_pms non-smoking ######################################
# no significant results from the ANCOM-BC data as there are NO hits

ancom.genus.old.pms = ancombc(phyloseq = genus.old.pms, # Raw counts
                              formula = 'age_disease_smoking', # The explanatory variable
                              p_adj_method = "fdr",
                              prv_cut=0.10, # Max proportion of zeros allowed per taxon
                              lib_cut = 1000, # Can filter out samples below minimum seq depth here
                              group = 'age_disease_smoking', # If you're including structural zeros below, you need this
                              struc_zero = T) # If true, any taxa present in only 1 category of 'group' will automatically be significant
head(ancom.genus.old.pms)
View(ancom.genus.old.pms)

# First we need to update the column names for each part of the results, because they're currently all the same.
colnames(ancom.genus.old.pms$res$lfc) = paste(colnames(ancom.genus.old.pms$res$lfc),'_beta',sep='')
colnames(ancom.genus.old.pms$res$se) = paste(colnames(ancom.genus.old.pms$res$se),'_se',sep='')
colnames(ancom.genus.old.pms$res$W) = paste(colnames(ancom.genus.old.pms$res$W),'_W',sep='')
colnames(ancom.genus.old.pms$res$p_val) = paste(colnames(ancom.genus.old.pms$res$p_val),'_p_val',sep='')
colnames(ancom.genus.old.pms$res$q_val) = paste(colnames(ancom.genus.old.pms$res$q_val),'_q_val',sep='')
colnames(ancom.genus.old.pms$res$diff_abn) = paste(colnames(ancom.genus.old.pms$res$diff_abn),'_diff_abn',sep='')
View(ancom.genus.old.pms)

#Family = age_disease_smoking
results = lapply(ancom.genus.old.pms$res,function(x) rownames_to_column(x,'age_disease_smoking')) %>% 
  lapply(as_tibble) %>% reduce(full_join)
results.sig = results %>% filter(`age_disease_smokingold, PMS, smoker_q_val_q_val`<0.05) # q_val = FDR-adjusted pval
hits = results.sig$age_disease_smoking 
hits = as.numeric(hits)
hits

# no significant results from the ANCOM-BC data as there are NO hits

# YPS = young_pms smoking v. young_pms non-smoking ##################################

ancom.genus.young.pms = ancombc(phyloseq = genus.young.pms, 
                                formula = 'age_disease_smoking', 
                                p_adj_method = "fdr",
                                prv_cut=0.10, 
                                lib_cut = 1000, 
                                group = 'age_disease_smoking', 
                                struc_zero = T) 
str(ancom.genus.young.pms)
View(ancom.genus.young.pms)

colnames(ancom.genus.young.pms$res$lfc) = paste(colnames(ancom.genus.young.pms$res$lfc),'_beta',sep='')
colnames(ancom.genus.young.pms$res$se) = paste(colnames(ancom.genus.young.pms$res$se),'_se',sep='')
colnames(ancom.genus.young.pms$res$W) = paste(colnames(ancom.genus.young.pms$res$W),'_W',sep='')
colnames(ancom.genus.young.pms$res$p_val) = paste(colnames(ancom.genus.young.pms$res$p_val),'_p_val',sep='')
colnames(ancom.genus.young.pms$res$q_val) = paste(colnames(ancom.genus.young.pms$res$q_val),'_q_val',sep='')
colnames(ancom.genus.young.pms$res$diff_abn) = paste(colnames(ancom.genus.young.pms$res$diff_abn),'_diff_abn',sep='')
View(ancom.genus.young.pms)

#Family = age_disease_smoking
results = lapply(ancom.genus.young.pms$res,function(x) rownames_to_column(x,'age_disease_smoking')) %>% 
  lapply(as_tibble) %>% reduce(full_join)
results.sig = results %>% filter(`age_disease_smokingyoung, PMS, smoker_q_val`<0.05) # q_val = FDR-adjusted pval
hits = results.sig$age_disease_smoking 
hits = as.numeric(hits)
hits

# WE GOT SIGNIFICANCE YAY!
# 2 20 21 36 41 44 52 59 61 66 67 68

genus_tss = genus.young.pms %>% 
                  transform('compositional')
selected_taxa <- taxa_names(genus.young.pms)[hits]
genus_tss <- prune_taxa(selected_taxa, genus_tss)
genus_tss_melt = genus_tss %>% 
                        psmelt()
bug = unique(genus_tss_melt$age_disease_smoking)
view(genus_tss_melt)

# why can I not see the age_disease_smoking variable

# This dataframe is super long - let's pivot it so that each microbe gets its own column.
# Note that sometimes you get an error saying that there is more than one value per row - this happens when bugs have the same family name, but belong to different taxonomic lineages. The best workaround is to add Order information to the Family column if this happens so that they can be differentiated. This commonly occurs when the family level isn't properly annotated by the software and is left as 'f__'.
family_tss_melt = genus_tss_melt %>% 
                    dplyr::select(-c(OTU,Domain:Family)) %>% # Remove all other levels; they'll interfere
                    pivot_wider(names_from = age_disease_smoking, values_from = Abundance)

# We'll use ggplot2 to make our plots. I don't want to write out 8 different plots and modify them individually, so I'm going to use a for loop here:
for(b in bug){
  #b = bug[1]
  
  #Just makes the column Fraction to uppercase
  p = family_tss_melt %>% 
    mutate(  = str_to_title( ))
  
  #generate random colors!
  colors = c(randomColor(1), randomColor(1))
  
  ggplot(p,aes(x = Fraction, y = p[[b]],fill=Fraction)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.2) +
    theme_classic(base_size = 16) +
    theme(legend.position='none') +
    scale_fill_manual(values = colors) +
    # Add the Padj (q) value from ANCOM-BC in Powerpoint, or hard code it in here
    xlab('Fraction') + ylab(paste(b,'(% Ab.)',sep=' '))
  print(p)
  # ggsave saves the last-generated plot
  ggsave(paste('Plots/tss_',b,'.jpeg',sep=''),height = 5,width = 5)
}

# young_healthy smoking v. young_healthy non-smoking

ancom.genus.young.healthy = ancombc(phyloseq = genus.young.healthy,
                                    formula = 'age_disease_smoking',
                                    p_adj_method = "fdr",
                                    prv_cut=0.10, 
                                    lib_cut = 1000, 
                                    group = 'age_disease_smoking', 
                                    struc_zero = T) 
str(ancom.genus.young.healthy)
View(ancom.genus.young.healthy)

# old_healthy smoking v. old_healthy non-smoking

ancom.genus.old.healthy = ancombc(phyloseq = genus.old.healthy, 
                                  formula = 'age_disease_smoking', 
                                  p_adj_method = "fdr",
                                  prv_cut=0.10, 
                                  lib_cut = 1000,
                                  group = 'age_disease_smoking', 
                                  struc_zero = T)
str(ancom.genus.old.healthy)
View(ancom.genus.old.healthy)

###############################################################################
# Step 4 - Visualizing ANCOM Results ##########################################









###############################################################################
# Step 5 - ALDeX2 #############################################################
###############################################################################

# ALDeX2: testing differential abundance across multiple conditions or groups and 
######### adjusts for within-sample variation.

library(ALDEx2)
set.seed(123) #for reproducibility

#################################################################################
## FAMILY LEVEL #################################################################
#################################################################################

# OPF = old_pms smoking v. old_pms non-smoking via family ##################################
s = family.old.pms@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = family.old.pms@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df_opf = aldex.glm(x) #old_pms_family
colnames(df_opf)
# "age_disease_smokingold, PMS, smoker:Est" = log_fold_change
# "age_disease_smokingold, PMS, smoker:pval" = p-value
# negative means that the ASVs are down-regulated

library(ggplot2)

log_fold_change <- df_opf$`age_disease_smokingold, PMS, smoker:Est`
pvals <- df_opf$`age_disease_smokingold, PMS, smoker:pval`
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
  labs(title = "Old_PMS, Smoking v. Non-Smoking ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# YPF = young_pms smoking v. young_pms non-smoking via family ##################################
s = family.young.pms@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = family.young.pms@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df_ypf = aldex.glm(x) # young_pms_family
colnames(df_ypf)

log_fold_change <- df_ypf$`age_disease_smokingyoung, PMS, smoker:Est`
pvals <- df_ypf$`age_disease_smokingyoung, PMS, smoker:pval`

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
  labs(title = "Young_PMS, Smoking v. Non-Smoking ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# YHF = young_healthy smoking v. young_healthy non-smoking via family ##########
s = family.young.healthy@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = family.young.healthy@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df_yhf = aldex.glm(x) # young_healthy_family
colnames(df_yhf)

log_fold_change <- df_yhf$`age_disease_smokingyoung, healthy, smoker:Est`
pvals <- df_yhf$`age_disease_smokingyoung, healthy, smoker:pval`

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
  labs(title = "Young_Healthy, Smoking v. Non-Smoking ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# NOTE. there seems to be TWO familiies of relevance

# OHF = old_healthy smoking v. old_healthy non-smoking via family ############
s = family.old.healthy@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = family.old.healthy@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df_ohf = aldex.glm(x) # old_healthy_family
colnames(df_ohf)

log_fold_change <- df_ohf$`age_disease_smokingold, healthy, smoker:Est`
pvals <- df_ohf$`age_disease_smokingold, healthy, smoker:pval`

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
  labs(title = "Old_Healthy, Smoking v. Non-Smoking ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# NOTE. there seems to be ONE family of relevance

################################################################################
## GENUS LEVEL #################################################################
################################################################################

# OPG = old_pms smoking v. old_pms non-smoking via genus ##################################
s = genus.old.pms@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = genus.old.pms@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df_opg = aldex.glm(x) #old_pms_genus
colnames(df_opg)

log_fold_change <- df_opg$`age_disease_smokingold, PMS, smoker:Est`
pvals <- df_opg$`age_disease_smokingold, PMS, smoker:pval`

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
  labs(title = "GENUS: Old_PMS, Smoking v. Non-Smoking ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# YPG = young_pms smoking v. young_pms non-smoking via genus ##################################
s = genus.young.pms@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = genus.young.pms@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df_ypg = aldex.glm(x) # young_pms_genus
colnames(df_ypg)

log_fold_change <- df_ypg$`age_disease_smokingyoung, PMS, smoker:Est`
pvals <- df_ypg$`age_disease_smokingyoung, PMS, smoker:pval`

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
  labs(title = "GENUS: Young_PMS, Smoking v. Non-Smoking ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# YHG = young_healthy smoking v. young_healthy non-smoking via genus ##########
s = genus.young.healthy@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = genus.young.healthy@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df_yhg = aldex.glm(x) # young_healthy_family
colnames(df_yhg)

log_fold_change <- df_yhg$`age_disease_smokingyoung, healthy, smoker:Est`
pvals <- df_yhg$`age_disease_smokingyoung, healthy, smoker:pval`

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
  labs(title = "GENUS: Young_Healthy, Smoking v. Non-Smoking ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")

# OHG = old_healthy smoking v. old_healthy non-smoking via genus ############
s = genus.old.healthy@sam_data %>% as.matrix() %>% as.data.frame() # metadata
m = model.matrix(~ age_disease_smoking, data = s) # converts to matrix 
o = genus.old.healthy@otu_table %>% as.matrix() %>% as.data.frame() # converts to dataframe
o = o[2:nrow(o),] %>% dplyr::select(all_of(rownames(m)))
x = aldex.clr(o,m,mc.samples=128)
df_ohg = aldex.glm(x) # old_healthy_genus
colnames(df_ohg)

log_fold_change <- df_ohg$`age_disease_smokingold, healthy, smoker:Est`
pvals <- df_ohg$`age_disease_smokingold, healthy, smoker:pval`

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
  labs(title = "GENUS: Old_Healthy, Smoking v. Non-Smoking ", 
       x = "Log2 Fold Change", 
       y = "-log10(p-value)") +
  theme_minimal() +
  theme(legend.position = "none")