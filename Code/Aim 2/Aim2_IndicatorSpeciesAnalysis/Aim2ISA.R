install.packages("indicspecies")

library(tidyverse)
library(phyloseq)
library(indicspecies)

load("MS_phyloseq.RData")


##group by age
## young indicates age is less than or equal to 55
## old indicates age is greater than 55
sample_data(MSphylo_filt)$agegroup <- ifelse(sample_data(MSphylo_filt)$age > 55, 'old', 'young')


##create column to group SPMS and PPMS for facet grid later on
sample_data(MSphylo_filt)$disease_group <- ifelse(
  sample_data(MSphylo_filt)$disease_course %in% c('PPMS', 'SPMS'), 'pMS', "Control"
)

##create new column for age+disease
sample_data(MSphylo_filt)$age_disease <- ifelse(
  sample_data(MSphylo_filt)$agegroup == "old" & sample_data(MSphylo_filt)$disease_group == "pMS", "old, PMS",
  ifelse(
    sample_data(MSphylo_filt)$agegroup == "young" & sample_data(MSphylo_filt)$disease_group == "pMS", "young, PMS",
    ifelse(
      sample_data(MSphylo_filt)$agegroup == "old" & sample_data(MSphylo_filt)$disease_group == 'Control', "old, healthy",
      ifelse(
        sample_data(MSphylo_filt)$agegroup == "young" & sample_data(MSphylo_filt)$disease_group == 'Control', "young, healthy",
        NA
      )
    )
  )
)

## create new column to add smoking status to age_disease
sample_data(MSphylo_filt)$age_disease_smoke <- ifelse(
  sample_data(MSphylo_filt)$age_disease == "old, PMS" & sample_data(MSphylo_filt)$smoke %in% c("smoker", "formersmoker"), "old, PMS, smoker",
  ifelse(
    sample_data(MSphylo_filt)$age_disease == "old, PMS" & sample_data(MSphylo_filt)$smoke == "nonsmoker", "old, PMS, nonsmoker",
    ifelse(
      sample_data(MSphylo_filt)$age_disease == "old, healthy" & sample_data(MSphylo_filt)$smoke %in% c("smoker", "formersmoker"), "old, healthy, smoker",
      ifelse(
        sample_data(MSphylo_filt)$age_disease == "old, healthy" & sample_data(MSphylo_filt)$smoke == "nonsmoker", "old, healthy, nonsmoker",
        ifelse(
          sample_data(MSphylo_filt)$age_disease == "young, PMS" & sample_data(MSphylo_filt)$smoke %in% c("smoker", "formersmoker"), "young, PMS, smoker",
          ifelse(
            sample_data(MSphylo_filt)$age_disease == "young, PMS" & sample_data(MSphylo_filt)$smoke == "nonsmoker", "young, PMS, nonsmoker",
            ifelse(
              sample_data(MSphylo_filt)$age_disease == "young, healthy" & sample_data(MSphylo_filt)$smoke %in% c("smoker", "formersmoker"), "young, healthy, smoker",
              ifelse(
                sample_data(MSphylo_filt)$age_disease == "young, healthy" & sample_data(MSphylo_filt)$smoke == "nonsmoker", "young, healthy, nonsmoker", NA
              )
            )
          )
        )
      )
    )
  )
)

     
get_variable(MSphylo_filt, c("age_disease_smoke"))

#### Indicator Species/Taxa Analysis ####
# glom to Genus
MS_genus <- tax_glom(MSphylo_filt, "Genus", NArm = FALSE)

#get relative abundance for indicator value
MS_genus_RA <- transform_sample_counts(MS_genus, fun=function(x) x/sum(x))

#ISA
isa_MS <- multipatt(t(otu_table(MS_genus_RA)), cluster = sample_data(MS_genus_RA)$`age_disease_smoke`)

taxtable <- tax_table(MSphylo_filt) %>% as.data.frame() %>% rownames_to_column(var="ASV")


# consider that your table is only going to be resolved up to the genus level, be wary of 
# anything beyond the glomed taxa level
isa_MS$sign %>%
  rownames_to_column(var="ASV") %>%
  left_join(taxtable) %>%
  filter(p.value<0.05) %>% View()