library(tidyverse)
library(cmapR)
library(org.Hs.eg.db)
library(rhdf5)

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

### Load gene expression data----

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


### Keep only signatures that have some duplicate sample----
duplicateSigs <- sigInfo %>% filter(dupl_counts>1)
duplicatesIndentity <- unique(duplicateSigs$duplIdentifier)

# Calculate GSEA distances

library(doRNG)

# run distances

# First get only data that have some duplicate signature.
# This will result in distance calculation between random
# samples and duplicates in a more computationaly traceable way
gex_dupls <- cmap[,unique(duplicateSigs$sig_id)]

# Specify the thresholds to check
# bottom and top regulated genes
thresholds <- c(30,50,100,200,300,400,
                500,600,700,800,900,1000)

# Initialize empty list for the results:
# Each element of the list (for each threshold)
# contains an NxN matrix with comparing all these
# samples. Each element of the matrix is the
# GSEA distance.
dist_all_dupls <- NULL

### SOS:
### RUN FIRST THE distance_scores.R
### SCRIPT TO LOAD THE FUNCTION!!!

### calculate distances: SEE distance_scores.R
# for information about the function inputs
dist_all_dupls <- foreach(thres = thresholds) %dorng% {
  distance_scores(num_table = gex_dupls ,
                  threshold_count = thres,names = colnames(gex_dupls))
}

# Transform list to array
distance <- do.call(cbind,dist_all_dupls)
distance <- array(distance,
                  c(dim=dim(dist_all_dupls[[1]]),length(dist_all_dupls)))

# Get the average distance across thresholds
mean_dist <- apply(distance, c(1,2), mean, na.rm = TRUE)
colnames(mean_dist) <- colnames(gex_dupls)
rownames(mean_dist) <- colnames(gex_dupls)

### Convert matrix into data frame
# Keep only unique (non-self) pairs
mean_dist[lower.tri(mean_dist,diag = T)] <- -100
dist <- reshape2::melt(mean_dist)
dist <- dist %>% filter(value != -100)

# Keep only specific columns for the original data frame
duplicateSigs <- duplicateSigs %>% 
  dplyr::select(sig_id,cmap_name,tas,nsample,duplIdentifier)

# Merge meta-data info and distances values
dist <- left_join(dist,duplicateSigs,by = c("Var1"="sig_id"))
dist <- left_join(dist,duplicateSigs,by = c("Var2"="sig_id"))
dist <- dist %>% mutate(is_duplicate = (duplIdentifier.x==duplIdentifier.y))

#Save results
saveRDS(dist,'../results/dupl_shRNA_mean_dist_q1replicates_5dupls_long.rds')

# Distance pre-processing on the pathways level----

# Perform pathway (gene set enrichment analysis)
keegEnrichResults <-  fastenrichment(sigInfo$sig_id,
                                     geneInfo$gene_id,
                                     cmap,
                                     enrichment_space = 'kegg',
                                     n_permutations=5000)
keggNES <- keegEnrichResults$NES$`NES KEGG`
keggpadj <- keegEnrichResults$Pval$`Pval KEGG`

saveRDS(keggNES,'../results/keggNES_shRNAs.rds')
saveRDS(keggpadj,'../results/keggpadj_shRNAs.rds')

#  
