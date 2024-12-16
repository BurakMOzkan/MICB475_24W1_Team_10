#!/usr/bin/env Rscript
# NOTE: to install phyloseq, please use the following code instead of the usual "install.packages" function:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
library(phyloseq)
library(ape) # importing trees
library(tidyverse)
library(vegan)
library(readr)

#### Load data ####
# Change file paths as necessary
MSmetafp <- "ms_metadata_filtered.tsv"
MSmeta <- read_delim(file=MSmetafp, delim="\t")

MSotufp <- "MS_feature_table.txt"
MSotu <- read_delim(file = MSotufp, delim="\t", skip=1)

MStaxfp <- "taxonomy.tsv"
MStax <- read_delim(MStaxfp, delim="\t")

MSphylotreefp <- "tree.nwk"
MSphylotree <- read.tree(MSphylotreefp)

#### Format OTU table ####
# OTU tables should be a matrix
# with rownames and colnames as OTUs and sampleIDs, respectively
# Note: tibbles do not allow rownames so if you imported with read_delim, change back

# save everything except first column (OTU ID) into a matrix
MSotu_mat <- as.matrix(MSotu[,-1])
# Make first column (#OTU ID) the rownames of the new matrix
rownames(MSotu_mat) <- MSotu$`#OTU ID`
# Use the "otu_table" function to make an OTU table
MSOTU <- otu_table(MSotu_mat, taxa_are_rows = TRUE) 
class(MSOTU)

#### Format sample metadata ####
# Save everything except sampleid as new data frame
MSsamp_df <- as.data.frame(MSmeta[,-1])
# Make sampleids the rownames
rownames(MSsamp_df)<- MSmeta$'sample-id'
# Make phyloseq sample data with sample_data() function
MSSAMP <- sample_data(MSsamp_df)
class(MSSAMP)

#### Formatting taxonomy ####
# Convert taxon strings to a table with separate taxa rank columns
MStax_mat <- MStax %>% select(-Confidence)%>%
  separate(col=Taxon, sep="; "
           , into = c("Domain","Phylum","Class","Order","Family","Genus","Species")) %>%
  as.matrix() # Saving as a matrix
# Save everything except feature IDs 
MStax_mat <- MStax_mat[,-1]
# Make sampleids the rownames
rownames(MStax_mat) <- MStax$`Feature ID`
# Make taxa table
MSTAX <- tax_table(MStax_mat)
class(MSTAX)

#### Create phyloseq object ####
# Merge all into a phyloseq object
MSPHYLOSEQ <- phyloseq(MSOTU, MSSAMP, MSTAX, MSphylotree)

MSphylo_filt <- subset_taxa(MSPHYLOSEQ, Domain == "d__Bacteria" & Class!="c__Chloroplast" & Family !="f__Mitochondria")


MSphylo_filt = subset_samples(MSphylo_filt, MSSAMP$`site_x`== 'San Francisco')

MSphylo_filt = subset_samples(MSphylo_filt, !is.na(smoke))

sample_data(MSphylo_filt)$age-group <- ifelse(sample_data(MSphylo_filt)$age > 55, 'old', 'young')

# Rarefy samples
# rngseed sets a random number. If you want to reproduce this exact analysis, you need
# to set rngseed the same number each time
# t transposes the table to use rarecurve function
# cex decreases font size
rarecurve(t(as.data.frame(otu_table(MSphylo_filt))), xlab='Sequencing Depth', ylab='Observed Features', cex=0.1)

MS_rare <- rarefy_even_depth(MSphylo_filt, rngseed = 1, sample.size = 6154)

save(MSphylo_filt, file="MS_phyloseq.RData")
save(MS_rare, file='MS_rare.RData')