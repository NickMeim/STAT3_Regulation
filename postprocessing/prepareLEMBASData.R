library(tidyverse)
library(cmapR)
library(org.Hs.eg.db)
library(rhdf5)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)


library(doFuture)

# parallel: set number of workers
cores <- 15
registerDoFuture()
plan(multiprocess,workers = cores)

### NEED TO CREATE A SIGNALING MODEL WITH INPUT THE GENE KNOCKOUTS NODES AND 
### FILTER SHRNAS NOT IN THE NETWORK
### FILTER CCLE GENES NOT IN THE NETWORK

### Load data and keep only well-inferred and landmark genes----
# Check L1000 documentation for information.
geneInfo <- read.delim('../data/geneinfo_beta.txt')
geneInfo <-  geneInfo %>% filter(feature_space != "inferred")
# Keep only protein-coding genes
geneInfo <- geneInfo %>% filter(gene_type=="protein-coding")

#Load signature info and split data to high quality replicates and low quality replicates
sigInfo <- read.delim('../data/siginfo_beta.txt')

# Create a proxy for quality of replicates
# Keep only samples with at least 3 replicates and that satisfy specific conditions.
# Check the LINCS2020 Release Metadata Field Definitions.xlsx file for 
# a complete description of each argument. It can be accessed online 
# or in the data folder.

sigInfo <- sigInfo %>% 
  mutate(quality_replicates = ifelse(is_exemplar_sig==1 & qc_pass==1 & nsample>=3,1,0))
# Filter only experiments of shrnas
sigInfo <- sigInfo %>% filter(pert_type=='trt_sh')
sigInfo <- sigInfo %>% filter(quality_replicates==1)

#sigInfo <- sigInfo %>% filter(tas>=0.25) # We keep all of them since they are statistically significant to have more hits

# Create identifier to signify duplicate
# signatures: meaning same drug, same dose,
# same time duration, same cell-type
sigInfo <- sigInfo %>% 
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))

sigInfo <- sigInfo %>% group_by(duplIdentifier) %>%
  mutate(dupl_counts = n()) %>% ungroup()

### Load CCLE data----
ccle <- t(data.table::fread('../data/CCLE/CCLE_expression.csv') %>% column_to_rownames('V1'))
ccle <- as.data.frame(ccle) %>% rownames_to_column('V1') %>% separate(V1,c('gene_id','useless'),sep=" ") %>%
  dplyr::select(-useless) %>% column_to_rownames('gene_id')
ccle <- as.data.frame(t(ccle)) %>% rownames_to_column('DepMap_ID')
sample_info <- data.table::fread('../data/CCLE/sample_info.csv') %>% dplyr::select(DepMap_ID,stripped_cell_line_name) %>%
  unique()
ccle <- left_join(ccle,sample_info) %>% dplyr::select(-DepMap_ID) %>%
  column_to_rownames('stripped_cell_line_name')
ccle <- ccle[which(rownames(ccle) %in% unique(sigInfo$cell_iname)),]
write.table(ccle, file = '../data/CCLE/preprocessed_ccle.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

### Create LEMBAS split for implementation with basal cell signal-----
sigInfo <- sigInfo %>% filter(cell_iname %in% rownames(ccle)) %>% unique()

### Split random 10fold validation----
library(caret)
total_samples <- sigInfo$sig_id
folds <- createFolds(total_samples, k = 10, list = TRUE, returnTrain = TRUE)

i <- 0
for (fold in folds){
  samples <- total_samples[fold]
  
  train_samples <- sigInfo %>% filter(sig_id %in% samples)
  val_samples <- sigInfo %>% filter(!(sig_id %in% samples))
  
  data.table::fwrite(as.data.frame(train_samples),paste0('../data/10fold_cross_validation/random/train_sample_',i,'.csv'),row.names = T)
  data.table::fwrite(as.data.frame(val_samples),paste0('../data/10fold_cross_validation/random/val_sample_',i,'.csv'),row.names = T)
  
  i <- i+1
}

### Split 3fold validation by keeping 4 cell-lines out every time----
library(caret)
total_samples <- unique(sigInfo$cell_iname)
folds <- createFolds(total_samples, k = 3, list = TRUE, returnTrain = TRUE)

i <- 0
for (fold in folds){
  samples <- total_samples[fold]
  
  train_samples <- sigInfo %>% filter(cell_iname %in% samples)
  val_samples <- sigInfo %>% filter(!(cell_iname %in% samples))
  
  data.table::fwrite(as.data.frame(train_samples),paste0('../data/3fold_cross_validation/cell_based/train_sample_',i,'.csv'),row.names = T)
  data.table::fwrite(as.data.frame(val_samples),paste0('../data/3fold_cross_validation/cell_based/val_sample_',i,'.csv'),row.names = T)
  
  i <- i+1
}
