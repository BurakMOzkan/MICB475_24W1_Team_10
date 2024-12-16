#this is the script for aim 1 alpha diversity analysis
#packages

library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(ggsignif)

#loading in rarefied data
load('MS_rare.Rdata')

##group by age
## young indicates age is less than or equal to 55
## old indicates age is greater than 55
sample_data(MS_rare)$agegroup <- ifelse(sample_data(MS_rare)$age > 55, 'old', 'young')


##create column to group SPMS and PPMS for facet grid later on
sample_data(MS_rare)$disease_group <- ifelse(
  sample_data(MS_rare)$disease_course %in% c('PPMS', 'SPMS'), 'pMS', "Control"
)

##create new column for age+disease
sample_data(MS_rare)$age_disease <- ifelse(
  sample_data(MS_rare)$agegroup == "old" & sample_data(MS_rare)$disease_group == "pMS", "older, pMS",
  ifelse(
    sample_data(MS_rare)$agegroup == "young" & sample_data(MS_rare)$disease_group == "pMS", "younger, pMS",
    ifelse(
      sample_data(MS_rare)$agegroup == "old" & sample_data(MS_rare)$disease_group == 'Control', "older, healthy",
      ifelse(
        sample_data(MS_rare)$agegroup == "young" & sample_data(MS_rare)$disease_group == 'Control', "younger, healthy",
        NA
      )
    )
  )
)

## now we have a new column which groups each sample based on their age group AND disease course
## this is what we will use for our analysis

## alpha diversity analysis

#plotting rarefied phyoseq object richness by 7 diff measures
#plot_richness(MS_rare, measures = c("Chao1", "Shannon", "InvSimpson"))

# if want to specify what metrics to use: 
# measures = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")

## plotting based on age and disease group with label as boxplot
# MS_gg_richness = plot_richness(MS_rare, measures = c("Chao1", "Shannon", "InvSimpson"), x = "age_disease") +
#   xlab("Age and Disease") + 
#   geom_boxplot()

# Add colours and decide as a group on colours for similar themes theme_classic to remove grey background,
# add layer of geom_points underboxplot line overlaid to see where points are,
# increase y (+ylim (0 lower limit)(upperlimit 0.5)) set dpi = 300 for ggsave

# MS_gg_richness_chao1 = plot_richness(MS_rare, measures = c("Chao1"), x = "agegroup") +
#   xlab("Age and Disease") + 
#   geom_boxplot() +
#   facet_grid(.~disease_group) +
#   theme
#   labs(title = "Chao1")
# 
# MS_gg_richness_shannon = plot_richness(MS_rare, measures = c("Shannon"), x = "agegroup") +
#   xlab("Age and Disease") + 
#   geom_boxplot() +
#   facet_grid(.~disease_group) +
#   labs(title = "Shannon")
# 
# MS_gg_richness_observed = plot_richness(MS_rare, measures = c("Observed"), x = "agegroup") +
#   xlab("Age and Disease") + 
#   geom_boxplot() +
#   facet_grid(.~disease_group) +
#   labs(title = "Observed")

MS_gg_richness_all = plot_richness(MS_rare, measures = c("Shannon", "Chao1", "Observed"), x = "age_disease") +
  xlab("Age and Disease") + 
  geom_boxplot() +
  geom_point() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
MS_gg_richness_all

MS_gg_richness_presentation = plot_richness(MS_rare, measures = c("Shannon"), x = "age_disease") +
  xlab("Age and Disease") + 
  geom_boxplot() +
  geom_point() +
  labs(title = "Aim 1 Richness") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#view the plots
MS_gg_richness_chao1
MS_gg_richness_shannon
MS_gg_richness_observed
MS_gg_richness_presentation
#save the plots
ggsave(filename = "MS_aim1_chao1.png",
       MS_gg_richness_chao1,
       height=4, width=10)

ggsave(filename = "MS_aim1_shannon.png",
       MS_gg_richness_shannon,
       height=4, width=10)

ggsave(filename = "MS_aim1_observed.png",
       MS_gg_richness_observed,
       height=4, width=10)

ggsave(filename = "MS_aim1_presentation.png",
       MS_gg_richness_presentation,
       height=4, width=10, dpi=300)

ggsave(filename = "MS_aim1_manu.png",
       MS_gg_richness_all,
       height=4, width=10, dpi=300)

## calculate faiths PD
phylo_dist <- pd(t(otu_table(MS_rare)), phy_tree(MS_rare),
                 include.root=F) 


# add PD to metadata table
sample_data(MS_rare)$PD <- phylo_dist$PD

# plot any age and disease group against the PD
MS_plot_pd <- ggplot(sample_data(MS_rare), aes(age_disease, PD)) + 
  geom_boxplot() +
  xlab("Age and Disease") +
  ylab("Phylogenetic Diversity")

MS_plot_pd

ggsave(filename = "MS_aim1_PD.png",
       MS_plot_pd,
       height=4, width=10)


### statistical tests for aim1 ###



alphadiv <- estimate_richness(MS_rare)
sam_dat <- sample_data(MS_rare)
samp_dat_wdiv <- data.frame(alphadiv, sam_dat)


##### SHANNON #####
##ANOVA LOG TRANS##
#set up linear model
lm_shannon_vs_age_disease <- lm(log(Shannon) ~ `age_disease`, dat = samp_dat_wdiv)
#calculate AOV
anova_shannon_age_disease <- aov(lm_shannon_vs_age_disease)
#summarize to determine if sig diffs
summary(anova_shannon_age_disease)
tukey_shannon <- TukeyHSD(anova_shannon_age_disease)
tukey_shannon

##KRUSKAL-WALLIS
kruskal_shannon <- kruskal.test(Shannon ~ `age_disease`, dat = samp_dat_wdiv)
kruskal_shannon
#chi-squared = 1.5825, df = 3, p-value = 0.6634


##### CHAO1 #####
##ANOVA LOG TRANS##
#set up linear model
lm_chao1_vs_age_disease <- lm(log(Chao1) ~ `age_disease`, dat = samp_dat_wdiv)
#calculate AOV
anova_chao1_age_disease <- aov(lm_chao1_vs_age_disease)
#summarize to determine if sig diffs
summary(anova_chao1_age_disease)
tukey_chao1 <- TukeyHSD(anova_chao1_age_disease)
tukey_chao1

##KRUSKAL-WALLIS
kruskal_chao1 <- kruskal.test(Chao1 ~ `age_disease`, dat = samp_dat_wdiv)
kruskal_chao1
#chi-sqaured = 5.804, df = 3, p-value = 0.1215


##### INVSIMPS #####
##ANOVA LOG TRANS##
#set up linear model
lm_invsimps_vs_age_disease <- lm(log(InvSimpson) ~ `age_disease`, dat = samp_dat_wdiv)
#calculate AOV
anova_invsimps_age_disease <- aov(lm_invsimps_vs_age_disease)
#summarize to determine if sig diffs
summary(anova_invsimps_age_disease)
tukey_invsimps <- TukeyHSD(anova_invsimps_age_disease)
tukey_invsimps

##KRUSKAL-WALLIS
kruskal_invsimps <- kruskal.test(InvSimpson ~ `age_disease`, dat = samp_dat_wdiv)
kruskal_invsimps
#chi-squared = 2.3196, df = 3, p-value = 0.5088


citation("ggVennDiagram")

