#packages
library(phyloseq)
library(ape)
library(tidyverse)
library(picante)
library(dplyr)
library(ggplot2)

##load MS data into environment
load("MS_rare.RData")

#check sample data
sample_data(MS_rare)

#check that only PPMS, SPMS, and control are selected
get_variable(MS_rare, c("disease_course"))

## Bin data based on age group
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


##Alpha diversity 

plot_richness(MS_rare) 

gg_richness <- plot_richness(MS_rare, x = "age_disease_smoking", measures = c("Observed", "Shannon","Chao1", "ACE", "Simpson", "InvSimpson", "Fisher")) +
  xlab("grouped individuals") +
  geom_boxplot()
gg_richness

ggsave(filename = "plot_richness.png"
       , gg_richness
       , height=10, width=20)


#################Further analysis for Shannon, Chao1, and Observed##########################
#shannon

gg_shannon <- plot_richness(MS_rare, x = "age_disease_smoking", measures = c("Shannon")) +
  xlab("grouped individuals") +
  geom_boxplot()
gg_shannon

ggsave(filename = "shannon.png"
       , gg_shannon
       , height=10, width=20)

alphadiversitymeasures <- estimate_richness(MS_rare)
alphadiversitymeasures

#add Shannon column into MS_rare data
sample_data(MS_rare)$Shannon <- alphadiversitymeasures$Shannon

get_variable(MS_rare, c("Shannon"))

#extract columns
samp_dat <- sample_data(MS_rare)
samp_dat_wdiv <- data.frame(samp_dat, alphadiversitymeasures)

##Statistical Analysis 

library(ggsignif)

####ANOVA
#set up linear model
lm_shannon <- lm(Shannon ~ `age_disease_smoking`, dat= samp_dat_wdiv)
# Calculate AOV
anova_shannon <- aov(lm_shannon)
# Summarize to determine if there are significant differences
summary(anova_shannon)
# Determine which groups are significant (tukey test for p value)
tukey_shannon <- TukeyHSD(anova_shannon)
tukey_shannon
#not significant, could be because this is a parametric test

####Non-parametric test: Kruskal-Wallis rank sum test
kruskal_shannon <- kruskal.test(Shannon ~ `age_disease_smoking`, data = samp_dat_wdiv)
kruskal_shannon
#p value is 0.08416 = not significant

# log the data and run ANOVA again to see if there is a difference
lm_log_shannon <- lm(log(Shannon) ~ `age_disease_smoking`, data=samp_dat_wdiv)
anova_log_shannon <- aov(lm_log_shannon)
summary(anova_log_shannon)
TukeyHSD(anova_log_shannon)
#no significance

#chao1
gg_chao1 <- plot_richness(MS_rare, x = "age_disease_smoking", measures = c("Observed")) +
  xlab("grouped individuals") +
  geom_boxplot()
gg_chao1

ggsave(filename = "chao1.png"
       , gg_chao1
       , height=10, width=20)

#add chao1 column into MS_rare data
sample_data(MS_rare)$Chao1 <- alphadiversitymeasures$Chao1

get_variable(MS_rare, c("Chao1"))

####ANOVA
#set up linear model
lm_Chao1 <- lm(Chao1 ~ `age_disease_smoking`, dat= samp_dat_wdiv)
# Calculate AOV
anova_Chao1 <- aov(lm_Chao1)
# Summarize to determine if there are significant differences
summary(anova_Chao1)
# Determine which groups are significant (tukey test for p value)
tukey_Chao1 <- TukeyHSD(anova_Chao1)
tukey_Chao1
#not significant

####Non-parametric test: Kruskal-Wallis rank sum test
kruskal_Chao1 <- kruskal.test(Chao1 ~ `age_disease_smoking`, data = samp_dat_wdiv)
kruskal_Chao1
#p value is 0.1199 = not significant
#data:  Chao1 by age_disease_smoking
#Kruskal-Wallis chi-squared = 11.456, df = 7, p-value = 0.1199

# log the data and run ANOVA again to see if there is a difference
lm_log_Chao1 <- lm(log(Chao1) ~ `age_disease_smoking`, data=samp_dat_wdiv)
anova_log_Chao1 <- aov(lm_log_Chao1)
summary(anova_log_Chao1)
TukeyHSD(anova_log_Chao1)
#no significant difference 

###observed with jitter plot
gg_observed <- plot_richness(MS_rare, x = "age_disease_smoking", measures = c("Observed")) +
  xlab("grouped individuals") +
  geom_boxplot() +
  geom_jitter(width= 0.1)
gg_observed

ggsave(filename = "observed.png"
       , gg_observed
       , height=10, width=20)

#shows similar to box but gives more specific values 
#jitter plot is easier to visualize a scatter plot
gg_observed1 <- plot_richness(MS_rare, x = "age_disease_smoking", measures = c("Observed")) +
  xlab("grouped individuals") +
  geom_jitter(width= 0.1)
gg_observed1

ggsave(filename = "jitter.observed.png"
       , gg_observed1
       , height=10, width=15)


#add observed column into MS_rare data
sample_data(MS_rare)$Observed <- alphadiversitymeasures$Observed

get_variable(MS_rare, c("Observed"))

######ANOVA
#set up linear model
lm_Observed <- lm(Observed ~ `age_disease_smoking`, dat= samp_dat_wdiv)
# Calculate AOV
anova_Observed <- aov(lm_Observed)
# Summarize to determine if there are significant differences
summary(anova_Observed)
# Determine which groups are significant (tukey test for p value)
tukey_Observed <- TukeyHSD(anova_Observed)
tukey_Observed
#not significant as expected LOL

####Non-parametric test: Kruskal-Wallis rank sum test
kruskal_Observed <- kruskal.test(Observed ~ `age_disease_smoking`, data = samp_dat_wdiv)
kruskal_Observed
#not significant
#Kruskal-Wallis chi-squared = 10.987, df = 7, p-value = 0.1392

# log the data and run ANOVA again to see if there is a difference
lm_log_Observed <- lm(log(Observed) ~ `age_disease_smoking`, data=samp_dat_wdiv)
anova_log_Observed <- aov(lm_log_Observed)
summary(anova_log_Observed)
TukeyHSD(anova_log_Observed)
#no significant difference between any group
