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

### Filter multiple time points
sigInfo <- sigInfo %>% mutate(timefree_id = paste0(cmap_name,'_',pert_idose,'_',cell_iname)) %>% group_by(timefree_id) %>%
  mutate(time_points=n_distinct(pert_itime)) %>% mutate(tas_max=max(tas)) %>%
  mutate(timekeep=ifelse(tas==tas_max | time_points==1,pert_itime,NA)) %>%  
  mutate(timekeep = unique(timekeep)[which(!is.na(unique(timekeep)))]) %>% ungroup() %>%
  filter(pert_itime==timekeep) %>% dplyr::select(-timekeep,-time_points,-timefree_id,-tas_max) %>% unique()

## Load omnipath
interactions <- read.delim('../preprocessing/preprocessed_data/FilteredOmnipath.tsv') %>% column_to_rownames('X')
## Filter shRNAs not in OmniPath
sigInfo <- sigInfo %>% filter(cmap_name %in% unique(c(interactions$source,interactions$target))) %>% unique()

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

### Infer transcription factors with Dorothea----
minNrOfGenes = 5
# Load requeired packages
library("dorothea")
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]

# Load GeX 
# Split sig_ids to run in parallel
sigIds <- unique(sigInfo$sig_id)
sigList <-  split(sigIds, 
                  ceiling(seq_along(sigIds)/ceiling(length(sigIds)/cores)))
# Parallelize parse_gctx function
# Path to raw data
ds_path <- '../../../../../L1000_2021_11_23/level5_beta_trt_sh_n238351x12328.gctx'
# rid is the gene entrez_id to find the gene in the data
# cid is the sig_id, meaning the sampe id
# path is the path to the data
parse_gctx_parallel <- function(path ,rid,cid){
  gctx_file <- parse_gctx(path ,rid = rid,cid = cid)
  return(gctx_file@mat)
}
# Parse the data file in parallel
cmap_gctx <- foreach(sigs = sigList) %dopar% {
  parse_gctx_parallel(ds_path ,
                      rid = unique(as.character(geneInfo$gene_id)),
                      cid = sigs)
}
cmap <-do.call(cbind,cmap_gctx)
df_annotation = data.frame(gene_id=rownames(cmap))
geneInfo$gene_id <- as.character(geneInfo$gene_id)
df_annotation <- left_join(df_annotation,geneInfo)
print(all(rownames(cmap)==df_annotation$gene_id))
rownames(cmap) <- df_annotation$gene_symbol

# Estimate TF activities
settings = list(verbose = TRUE, minsize = minNrOfGenes)
TF_activities = run_viper(cmap, dorotheaData, options =  settings)

TF_activities <- 1/(1+exp(-TF_activities))
TF_activities <- t(TF_activities)
hist(TF_activities)
TF_activities <- TF_activities[,which(colnames(TF_activities) %in% unique(c(interactions$source,interactions$target)))]
write.table(TF_activities, file = '../results/filtered_shrnas_tf_activities.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)


### Create signaling model------------------

### Targeted TFS
tfs_tageted <- colnames(TF_activities)[which(colnames(TF_activities) %in% unique(sigInfo$cmap_name))]
tfs_tageted <- as.data.frame(tfs_tageted)
colnames(tfs_tageted) <- 'Entry'
write.table(tfs_tageted, file = '../preprocessing/preprocessed_data/targetd_tfs.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

### FILTER SHRNAS AND TFS NOT IN THE TRIMMED NETWORK
pkn <- read.delim('../preprocessing/preprocessed_data/PKN-Model.tsv')
sigInfo <- sigInfo %>% filter(cmap_name %in% unique(c(pkn$source,pkn$target)))
TF_activities <- TF_activities[,which(colnames(TF_activities) %in% unique(c(pkn$source,pkn$target)))]
write.table(TF_activities, file = '../results/trimmed_shrnas_tf_activities.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

### FILTER CCLE GENES NOT IN THE NETWORK ---> those not in the trimmed omnipath.
ccle <- ccle[,which(colnames(ccle) %in% unique(c(pkn$source,pkn$target)))]
ccle <- ccle[which(rownames(ccle) %in% unique(sigInfo$cell_iname)),]
write.table(ccle, file = '../data/CCLE/trimmed_ccle.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

### Split random 10fold validation----
sigInfo <- sigInfo %>% filter(cell_iname %in% rownames(ccle)) %>% unique()

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
  
  print(paste0('Train cells:',length(unique(train_samples$cell_iname))))
  print(paste0('Val cells:',length(unique(val_samples$cell_iname))))
  
  data.table::fwrite(as.data.frame(train_samples),paste0('../data/3fold_cross_validation/cell_based/train_sample_',i,'.csv'),row.names = T)
  data.table::fwrite(as.data.frame(val_samples),paste0('../data/3fold_cross_validation/cell_based/val_sample_',i,'.csv'),row.names = T)
  
  i <- i+1
}

### Knock-outs are simulated as putting -5 in the input node----
# Create a general condition input data (there are not multiple doses (only NAs and a specific dose))
# Probably finding the cmap_name while training and saying -5 there is sufficient
cellInfo <- sigInfo %>% dplyr::select(sig_id,cell_iname) %>% unique()
sigInfo <- sigInfo %>% dplyr::select(sig_id,cmap_name) %>% unique()
sigInfo <- sigInfo %>% mutate(value=-5) %>% spread('cmap_name','value')
sigInfo[is.na(sigInfo)] <- 0
write.table(sigInfo, file = '../preprocessing/preprocessed_data/all_filtered_Kds.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
cellInfo <- cellInfo %>% mutate(value=1) %>% spread('cell_iname','value')
cellInfo[is.na(cellInfo)] <- 0
write.table(cellInfo, file = '../preprocessing/preprocessed_data/all_filtered_cells.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
