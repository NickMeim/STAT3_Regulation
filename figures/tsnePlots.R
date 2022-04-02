# Plot STAT3 in t-SNE space and then 1+9 random duplicates 
# Visualize the latent space in 2D
library(tidyverse)
library(ggplot2)
library(ggpubr)

#### Load meta-data of the L1000 platform
### It is the same as in pre-procesing

# Check L1000 documentation for information.

#Load signature info and split data to 
#high quality replicates and low quality replicates
sigInfo <- read.delim('../data/siginfo_beta.txt')

# Create a proxy for quality of replicates
# Keep only samples with at least 3 replicates and that satisfy specific conditions.
# Check the LINCS2020 Release Metadata Field Definitions.xlsx file for 
# a complete description of each argument. It can be accessed online 
# or in the data folder.

sigInfo <- sigInfo %>% 
  mutate(quality_replicates = ifelse(is_exemplar_sig==1 & qc_pass==1 & nsample>=3,1,0))
# Filter only experiments of drugs
sigInfo <- sigInfo %>% filter(pert_type=='trt_sh')
sigInfo <- sigInfo %>% filter(quality_replicates==1)

# Create identifier to signify duplicate
# signatures: meaning same drug, same dose,
# same time duration, same cell-type
sigInfo <- sigInfo %>%
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))

sigInfo <- sigInfo %>% group_by(duplIdentifier) %>% 
  mutate(dupl_counts = n()) %>% ungroup()

# Load embeddings
emb <- read.csv('../../shrna_embs1024.csv')
emb <- emb %>% column_to_rownames('X')
emb <- as.matrix(emb)
gc()

# Visualize t-SNE STAT3 vs others


# Keep some duplicates to be computationally feasible
duplicateSigs <- sigInfo %>% filter(dupl_counts>=5)
duplicatesIndentity <- unique(duplicateSigs$duplIdentifier)

emb <- emb[as.character(unique(duplicateSigs$sig_id)),]
gc()

