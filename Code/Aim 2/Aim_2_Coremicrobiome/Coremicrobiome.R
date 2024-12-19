library(phyloseq)
library(microbiome)
library(ggVennDiagram)
library(tibble)

#load the phyloseq object
phylo <- MSphylo_filt
phyloseq_object <- transform_sample_counts(phylo, fun=function(x) x/sum(x))

#Subset samples based on age, smoking, and disease course
young_nonsmokers_control <- subset_samples(phyloseq_object, age <= 55 & disease_course == "Control" & smoke == "nonsmoker")
old_nonsmokers_control<- subset_samples(phyloseq_object, age >= 56 & disease_course == "Control" & smoke == "nonsmoker")

old_nonsmokers_ms<- subset_samples(phyloseq_object, age >= 56 & disease == "MS" & smoke == "nonsmoker")
young_nonsmokers_ms <- subset_samples(phyloseq_object, age <= 55 & disease == "MS" & smoke == "nonsmoker")

young_smokers_control <- subset_samples(phyloseq_object, age <= 55 & disease_course == "Control" & (smoke == "formersmoker" | smoke == "smoker"))
old_smokers_control<- subset_samples(phyloseq_object, age >= 56 & disease_course == "Control" & (smoke == "formersmoker" | smoke == "smoker"))

old_smokers_ms<- subset_samples(phyloseq_object, age >= 56 & disease == "MS" & (smoke == "formersmoker" | smoke == "smoker"))
young_smokers_ms <- subset_samples(phyloseq_object, age <= 55 & disease == "MS" & (smoke == "formersmoker" | smoke == "smoker"))

#Conduct coremembers analusis, set detection to 0.001 and prevalence to 0.5
young_nonsmokers_control_core <- core_members(young_nonsmokers_control, detection=0.001, prevalence=0.5)
old_nonsmokers_control_core <- core_members(old_nonsmokers_control, detection=0.001, prevalence=0.5)
old_nonsmokers_ms_core <- core_members(old_nonsmokers_ms, detection=0.001, prevalence=0.5)
young_nonsmokers_ms_core <- core_members(young_nonsmokers_ms, detection=0.001, prevalence=0.5)
young_smokers_control_core <- core_members(young_smokers_control, detection=0.001, prevalence=0.5)
old_smokers_control_core <- core_members(old_smokers_control, detection=0.001, prevalence=0.5)
old_smokers_ms_core <- core_members(old_smokers_ms, detection=0.001, prevalence=0.5)
young_smokers_ms_core <- core_members(young_smokers_ms, detection=0.001, prevalence=0.5)

#Produce venn diagrams 
young_venn <- ggVennDiagram(x=list(young_smokers_control_core, young_smokers_ms_core, young_nonsmokers_control_core, young_nonsmokers_ms_core), 
                            category.names = c("Control Smokers", "pMS Smokers", "Control Nonsmokers", "pMS Nonsmokers")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  ggtitle("Core Microbiomes of Younger Cohort") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) # Center and style the title )
old_venn <- ggVennDiagram(x=list(old_smokers_control_core, old_smokers_ms_core, old_nonsmokers_control_core, old_nonsmokers_ms_core), 
                          category.names = c("Control Smokers", "pMS Smokers", "Control Nonsmokers", "pMS Nonsmokers")) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  ggtitle("Core Microbiomes of Older Cohort") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) # Center and style the title )


young_venn
old_venn

