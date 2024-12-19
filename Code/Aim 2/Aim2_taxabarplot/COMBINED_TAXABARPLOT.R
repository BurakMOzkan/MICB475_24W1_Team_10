library(phyloseq)
library(indicspecies)
library(microbiome)
library(tidyverse)

phylo <- MSphylo_filt
physeq_phylum <- tax_glom(phylo, taxrank = "Class")
phyloseq_object <- transform_sample_counts(physeq_phylum, fun=function(x) x/sum(x))

#bin by age
phyloseq_object@sam_data$agegroup <- cut(phyloseq_object@sam_data$age, 
                                         breaks = c(25, 55, Inf), 
                                         labels = c("young", "old"),
                                         right = FALSE)

#bin by age and disease 
sample_data(phyloseq_object)$age_disease <- ifelse(
  sample_data(phyloseq_object)$agegroup == "old" & sample_data(phyloseq_object)$disease_course %in% c('PPMS', 'SPMS'), "old, PMS",
  ifelse(
    sample_data(phyloseq_object)$agegroup == "young" & sample_data(phyloseq_object)$disease_course %in% c('PPMS', 'SPMS'), "young, PMS",
    ifelse(
      sample_data(phyloseq_object)$agegroup == "old" & sample_data(phyloseq_object)$disease_course == 'Control', "old, healthy",
      ifelse(
        sample_data(phyloseq_object)$agegroup == "young" & sample_data(phyloseq_object)$disease_course == 'Control', "young, healthy",
        NA ))))

#bin by age, disease, and smoking 
sample_data(phyloseq_object)$age_disease_smoking <- ifelse(
  sample_data(phyloseq_object)$age_disease == "old, healthy" & sample_data(phyloseq_object)$smoke %in% c('formersmoker', 'smoker'), "older, healthy, smoker", 
  ifelse(
    sample_data(phyloseq_object)$age_disease == "young, healthy" & sample_data(phyloseq_object)$smoke %in% c('formersmoker', 'smoker'), "younger, healthy, smoker",
    ifelse(
      sample_data(phyloseq_object)$age_disease == "old, healthy" & sample_data(phyloseq_object)$smoke == 'nonsmoker', "older, healthy, nonsmoker",
      ifelse(
        sample_data(phyloseq_object)$age_disease == "young, healthy" & sample_data(phyloseq_object)$smoke == 'nonsmoker', "younger, healthy, nonsmoker",
        ifelse(
          sample_data(phyloseq_object)$age_disease == "old, PMS" & sample_data(phyloseq_object)$smoke == 'nonsmoker', "older, pMS, nonsmoker",
          ifelse(
            sample_data(phyloseq_object)$age_disease == "young, PMS" & sample_data(phyloseq_object)$smoke == 'nonsmoker', "younger, pMS, nonsmoker",
            ifelse(
              sample_data(phyloseq_object)$age_disease == "old, PMS" & sample_data(phyloseq_object)$smoke %in% c('formersmoker', 'smoker'), "older, pMS, smoker",
              ifelse(
                sample_data(phyloseq_object)$age_disease == "young, PMS" & sample_data(phyloseq_object)$smoke %in% c('formersmoker', 'smoker'), "younger, pMS, smoker", NA ))))))))

physeq_phylum_combined <- merge_samples(phyloseq_object, group = "age_disease_smoking")

physeq_genus_melt <- psmelt(physeq_phylum_combined)

# 5. Create the Stacked Bar Plot
 ggplot(physeq_genus_melt, aes(x = Sample, y = Abundance, fill = Class)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(y = "Relative Abundance", fill = "Class") +
  scale_y_continuous(labels = scales::percent) +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 16))

ggsave("combined_taxabarplot", plot = save_plot, device = "png",  width = 10, height = 7.5, units = "in")
