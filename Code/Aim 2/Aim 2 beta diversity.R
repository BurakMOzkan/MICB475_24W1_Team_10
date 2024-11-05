#!/usr/bin/env Rscript
library(tidyverse)
library(phyloseq)
library(vegan)

##load MS data into environment
load("MS_rare.RData")

#check sample data
sample_data(MS_rare)

#check that only PPMS, SPMS, and control are selected
get_variable(MS_rare, c("disease_course", "smoke"))

## bin data based on age group
MS_rare@sam_data$agegroup <- cut(MS_rare@sam_data$age, 
                                 breaks = c(25, 55, Inf), 
                                 labels = c("young", "old"),
                                 right = FALSE)
MS_rare@sam_data$agegroup

## Make a new column for the chosen variables 
##create new column for age and disease
sample_data(MS_rare)$age_disease <- ifelse(
  sample_data(MS_rare)$agegroup == "old" & sample_data(MS_rare)$disease_course %in% c('PPMS', 'SPMS'), "old, PMS",
  ifelse(
    sample_data(MS_rare)$agegroup == "young" & sample_data(MS_rare)$disease_course %in% c('PPMS', 'SPMS'), "young, PMS",
    ifelse(
      sample_data(MS_rare)$agegroup == "old" & sample_data(MS_rare)$disease_course == 'Control', "old, healthy",
      ifelse(
        sample_data(MS_rare)$agegroup == "young" & sample_data(MS_rare)$disease_course == 'Control', "young, healthy",
        NA ))))
get_variable(MS_rare, c("age_disease"))

##Make a column for age, disease, and smoking status
sample_data(MS_rare)$age_disease_smoking <- ifelse(
  sample_data(MS_rare)$age_disease == "old, healthy" & sample_data(MS_rare)$smoke %in% c('formersmoker', 'smoker'), "old, healthy, smoker",
  ifelse(
    sample_data(MS_rare)$age_disease == "young, healthy" & sample_data(MS_rare)$smoke %in% c('formersmoker', 'smoker'), "young, healthy, smoker",
    ifelse(
      sample_data(MS_rare)$age_disease == "old, healthy" & sample_data(MS_rare)$smoke == 'nonsmoker', "old, healthy, nonsmoker",
      ifelse(
        sample_data(MS_rare)$age_disease == "young, healthy" & sample_data(MS_rare)$smoke == 'nonsmoker', "young, healthy, nonsmoker",
        ifelse(
          sample_data(MS_rare)$age_disease == "old, PMS" & sample_data(MS_rare)$smoke == 'nonsmoker', "old, PMS, nonsmoker",
          ifelse(
            sample_data(MS_rare)$age_disease == "young, PMS" & sample_data(MS_rare)$smoke == 'nonsmoker', "young, PMS, nonsmoker",
            ifelse(
              sample_data(MS_rare)$age_disease == "old, PMS" & sample_data(MS_rare)$smoke %in% c('formersmoker', 'smoker'), "old, PMS, smoker",
              ifelse(
                sample_data(MS_rare)$age_disease == "young, PMS" & sample_data(MS_rare)$smoke %in% c('formersmoker', 'smoker'), "young, PMS, smoker", NA ))))))))
get_variable(MS_rare, c("age_disease_smoking"))

samp_dat_wdiv <- data.frame(sample_data(MS_rare), estimate_richness(MS_rare))


#Metric: unifrac
dm_unifrac <- UniFrac(MS_rare, weighted=TRUE)
# plot the above as an ordination to a PCoA plot
ord.unifrac <- ordinate(MS_rare, method="PCoA", distance="unifrac")
plot_ordination(MS_rare, ord.unifrac, color="age_disease_smoking")+
  stat_ellipse(type = "norm")
  labs(title="Aim 2 Beta Diversity, Unweighted Unifrac")
#save file
  ggsave(filename = "plot_uw_unifrac.png"
         , gg_richness
         , height=10, width=20)




# Metric: PERMANOVA weighted unifrac
adonis2(dm_unifrac ~ age_disease_smoking, data=samp_dat_wdiv)
#plot
plot_ordination(MS_rare, ord.unifrac, color="age_disease_smoking")+
  stat_ellipse(type = "norm")+
  labs(title="Aim 2 Beta Diversity, Weighted  Unifrac")
#save file
ggsave(filename = "plot_w_unifrac.png"
       , gg_richness
       , height=10, width=20)

  
  
# Metric: Bray
dm_bray <- vegdist(t(otu_table(MS_rare)), method="bray")
adonis2(dm_bray ~ age_disease_smoking, data=samp_dat_wdiv)
#plot
plot_ordination(MS_rare, dm_bray, color='age_disease_smoking')+
  stat_ellipse(type = "norm")
  labs(title="Aim 2 Beta Diversity, Bray")
#save file
  ggsave(filename = "plot_bray.png"
         , gg_richness
         , height=10, width=20)




#Matric Jaccard
dm_jaccard <- vegdist(t(otu_table(MS_rare)), method="jaccard")
adonis2(dm_jaccard ~ age_disease_smoking, data=samp_dat_wdiv)
#plot
plot_ordination(MS_rare, dm_jaccard, color='age_disease_smoking')+
  stat_ellipse(type = "norm")
  labs(title="Aim 2 Beta Diversity, Jaccard")
#save file
  ggsave(filename = "plot_jaccard.png"
         , gg_richness
         , height=10, width=20)

#statistical analysis
anova(lm(disease_course ~ smoke, data=MS_rare))








