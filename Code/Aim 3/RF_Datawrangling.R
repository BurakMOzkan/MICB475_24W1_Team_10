library(phyloseq)
library(microbiome)
library(tidyverse)
library(readr)
library(dplyr)

phylo <- MSphylo_filt


family_of_interest <- c("f__Akkermansiaceae", "f__Peptostreptococcales-Tissierellales", 
                        "f__Butyricicoccaceae", "f__Oscillospiraceae")

genus_of_interest <- c("g__Oscillibacter", "g__Akkermansia", "g__Faecalibacterium", 
                       "g__Erysipelotrichaceae_UCG", "g__Roseburia", "g__Christensenellaceae_R-7_group", 
                       "g__Blautia", "g__Fusicatenibacter", "g__Agathobacter", "g__Frisingicoccus")

metadata_table <- phylo@sam_data

ps_for_rf_genus = phylo %>% 
  microbiome::transform('clr') %>%  # Works on phyloseq objects
  subset_taxa(Genus %in% genus_of_interest) %>%
  tax_glom("Family") %>% psmelt() %>%
  select(-OTU, -c(Domain:Order)) %>% 
  pivot_wider(names_from = Family, values_from = Abundance) %>%
  select(Sample,age, smoke, disease, contains('f__'))

ps_for_rf_family = phylo %>% 
  microbiome::transform('clr') %>%  # Works on phyloseq objects
  subset_taxa(Family %in% all_families_final) %>%
  tax_glom("Family") %>% psmelt() %>%
  select(-OTU, -c(Domain:Order)) %>% 
  pivot_wider(names_from = Family, values_from = Abundance) %>%
  select(Sample, contains('f__'))

combined <- left_join(ps_for_rf_genus, ps_for_rf_family, by = "Sample") 

ps_for_rf_genus <- sapply(ps_for_rf_genus$list_column, `[`, 1)

ps_for_rf_genus <- as.data.frame(ps_for_rf_genus)

write.csv(ps_for_rf_family, file = "MS_dataforRF_allgenusallage.csv", row.names = FALSE)
