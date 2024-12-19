library(phyloseq)
library(microbiome)
library(tidyverse)
library(readr)
library(dplyr)

phylo <- MSphylo_filt

#Determine list of family/genus we want to incorporate into 
family_of_interest <- c("f__Akkermansiaceae", "f__Peptostreptococcales-Tissierellales", 
                        "f__Butyricicoccaceae", "f__Oscillospiraceae")

genus_of_interest <- c("g__Oscillibacter", "g__Akkermansia", "g__Faecalibacterium", 
                       "g__Erysipelotrichaceae_UCG", "g__Roseburia", "g__Christensenellaceae_R-7_group", 
                       "g__Blautia", "g__Fusicatenibacter", "g__Agathobacter", "g__Frisingicoccus")

#Create table of metadata
metadata_table <- phylo@sam_data

#Create dataframe of genus/family values
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

#Combine family and genus level data together
combined <- left_join(ps_for_rf_genus, ps_for_rf_family, by = "Sample") 

#Save File 
write.csv(ps_for_rf_family, file = "MS_dataforRF_allgenusallage.csv", row.names = FALSE)
