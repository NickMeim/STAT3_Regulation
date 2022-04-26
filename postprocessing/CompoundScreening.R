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


### Therapeutic screening: Load data----

# Check L1000 documentation for information.
geneInfo <- read.delim('../data/geneinfo_beta.txt')
geneInfo <-  geneInfo %>% filter(feature_space != "inferred")
# Keep only protein-coding genes
geneInfo <- geneInfo %>% filter(gene_type=="protein-coding")

# Screen drugs based on the GSEA distance with the STAT3 knockdown
#Load signature info and split data to high quality replicates and low quality replicates
sigInfo <- read.delim('../data/siginfo_beta.txt')

# Create a proxy for quality of replicates
# Keep only samples with at least 3 replicates and that satisfy specific conditions.
# Check the LINCS2020 Release Metadata Field Definitions.xlsx file for 
# a complete description of each argument. It can be accessed online 
# or in the data folder.

sigInfo <- sigInfo %>% 
  mutate(quality_replicates = ifelse(is_exemplar_sig==1 & qc_pass==1 & nsample>=3,1,0))

# Filter STAT3 shrnas
shRNA <- sigInfo %>% filter(pert_type=='trt_sh')
shRNA <- shRNA %>% filter(quality_replicates==1)
shRNA <-  shRNA %>% filter(cmap_name=='STAT3' | tas>=0.3)
# Keep only relevant candidates
nodeAttr <- read.delim('../results/releventNodeAttr.txt')
knockdowns <- unique(nodeAttr$node)
shRNA <- shRNA %>% filter(cmap_name %in% knockdowns)

# Filter drugs
sigInfo <- sigInfo %>% filter(pert_type=='trt_cp')
sigInfo <- sigInfo %>% filter(quality_replicates==1)
sigInfo <- sigInfo %>% filter(tas>=0.3)

# Create identifier to signify duplicate
# signatures: meaning same drug, same dose,
# same time duration, same cell-type
sigInfo <- sigInfo %>% 
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))
sigInfo <- sigInfo %>% group_by(duplIdentifier) %>%
  mutate(dupl_counts = n()) %>% ungroup()
shRNA <- shRNA %>% 
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))
shRNA <- shRNA %>% group_by(duplIdentifier) %>%
  mutate(dupl_counts = n()) %>% ungroup()

### Therapeutic screening: Load gene expression data----

# rid is the gene entrez_id to find the gene in the data
# cid is the sig_id, meaning the sampe id
# path is the path to the data
parse_gctx_parallel <- function(path ,rid,cid){
  gctx_file <- parse_gctx(path ,rid = rid,cid = cid)
  return(gctx_file@mat)
}
# Load STAT3 and shRNA candidates GeX to use it later
# in distance calculations
# Path to raw shRNA data
shRNA_path <- '../../../../../L1000_2021_11_23/level5_beta_trt_sh_n238351x12328.gctx'
cmap_shRNA <- parse_gctx(shRNA_path ,
                        rid = unique(as.character(geneInfo$gene_id)),
                        cid = unique(as.character(shRNA$sig_id)))
cmap_shRNA <- cmap_shRNA@mat

# Now read GeX for drugs
# Split sig_ids to run in parallel
sigIds <- unique(sigInfo$sig_id)
sigList <-  split(sigIds, 
                  ceiling(seq_along(sigIds)/ceiling(length(sigIds)/cores)))

# Parallelize parse_gctx function

# Path to raw data
ds_path <- '../../../../../L1000_2021_11_23/level5_beta_trt_cp_n720216x12328.gctx'

# Parse the data file in parallel
cmap_gctx <- foreach(sigs = sigList) %dopar% {
  parse_gctx_parallel(ds_path ,
                      rid = unique(as.character(geneInfo$gene_id)),
                      cid = sigs)
}
cmap <-do.call(cbind,cmap_gctx)

# Bind together GeX of drugs and STAT3 shRNA (together with also candidates shRNA) 
# to calculate all pairwise distances
cmap_combined <- cbind(cmap,cmap_shRNA)

### Therapeutic screening: Perform GSEA to calculate pathways enrichment----
sigInfo_combined <- rbind(sigInfo,shRNA)
keegEnrichResults <-  fastenrichment(sigInfo_combined$sig_id,
                                     geneInfo$gene_id,
                                     cmap_combined,
                                     enrichment_space = 'kegg',
                                     n_permutations=5000,
                                     pval_adjustment=F)

keggNES <- keegEnrichResults$NES$`NES KEGG`
keggpval <- keegEnrichResults$Pval$`Pval KEGG`

saveRDS(keggNES,'../results/keggNES_combined_drugs_shrna.rds')
saveRDS(keggpval,'../results/keggpval_combined_drugs_shrna.rds')


# Calculate distance of candidates from STAT3
keggNES <- readRDS('../results/keggNES_combined_drugs_shrna.rds')
keggNES <-  keggNES[,as.character(sigInfo$sig_id)]

# Calculate gsea distances 
library(doRNG)
# run distances
thresholds <- c(5,10,20,30,40,50)
dist_all_dupls <- NULL
print('Begin calculating GSEA distance...')
### calculate distances
dist_all_dupls <- foreach(thres = thresholds) %dorng% {
  distance_scores(num_table = keggNES ,threshold_count = thres,names = colnames(keggNES))
}
distance <- do.call(cbind,dist_all_dupls)
distance <- array(distance,c(dim=dim(dist_all_dupls[[1]]),length(dist_all_dupls)))
mean_dist <- apply(distance, c(1,2), mean, na.rm = TRUE)
colnames(mean_dist) <- colnames(keggNES)
rownames(mean_dist) <- colnames(keggNES)
print('Begin saving GSEA distance...')
#saveRDS(mean_dist,'../results/cmap_mean_dist_kegg_shrna_statcandidates_tasfiltered.rds')

# Keep good candidates
mean_dist <- readRDS('../results/cmap_mean_dist_kegg_shrna_statcandidates_tasfiltered.rds')
gc()

### Convert matrix into data frame
# Keep only unique (non-self) pairs
mean_dist[lower.tri(mean_dist,diag = T)] <- -100
dist <- reshape2::melt(mean_dist)
dist <- dist %>% filter(value != -100)

# Merge meta-data info and distances values
dist <- left_join(dist,sigInfo %>% dplyr::select(sig_id,cmap_name,cell_iname,duplIdentifier),by = c("Var1"="sig_id"))
dist <- left_join(dist,sigInfo %>% dplyr::select(sig_id,cmap_name,cell_iname,duplIdentifier),by = c("Var2"="sig_id"))
dist <- dist %>% mutate(is_duplicate = (duplIdentifier.x==duplIdentifier.y))
dist$value <- dist$value/2

# Keep only STAT3 with other candidate pairs
dist <- dist %>% filter(cmap_name.x=='STAT3' | cmap_name.y=='STAT3')
dist <- dist %>% filter(!(cmap_name.x=='STAT3' & cmap_name.y=='STAT3'))
dist <- dist %>% filter(cell_iname.x==cell_iname.y)
dist <- dist %>% mutate(cell_iname = cell_iname.x) %>% dplyr::select(-cell_iname.x,-cell_iname.y) %>% unique()

# First aggregate duplicates
dist <- dist %>% mutate(pairID=ifelse(cmap_name.x!='STAT3',duplIdentifier.x,duplIdentifier.y)) %>%
  group_by(pairID) %>% mutate(med_value=median(value)) %>% ungroup()
dist$value <- dist$med_value
dist <- dist %>% dplyr::select(-med_value,-Var1,-Var2,-duplIdentifier.x,-duplIdentifier.y) %>% unique()
#df <- dist %>% group_by(pairID) %>% summarise(n())

# Get neighbors per cell
cells <- unique(dist$cell_iname)

#dist <- dist %>% mutate(similar=ifelse(value<=0.3,1,0))
dist_filtered <- dist %>% filter(value<=0.3) # found from duplicates distribution

# Get counts of how many times each candidate was found as a neighbor
genes <- unique(c(dist_filtered$cmap_name.x,dist_filtered$cmap_name.y))
genes <- genes[-which(genes=='STAT3')]

dist_filtered <- dist_filtered %>% mutate(proportion=0)
for (i in 1:length(genes)){
  inds <- which(dist_filtered$cmap_name.x==genes[i] | dist_filtered$cmap_name.y==genes[i])
  df <- dist %>% filter(cmap_name.x==genes[i] | cmap_name.y==genes[i]) %>% dplyr::select(-cmap_name.x,-cmap_name.y) %>% unique()
  df_filtered <- dist_filtered %>% filter(cmap_name.x==genes[i] | cmap_name.y==genes[i]) %>% dplyr::select(-cmap_name.x,-cmap_name.y) %>% unique()
  proportion <- length(unique(df_filtered$cell_iname))/length(unique(df$cell_iname))
  dist_filtered$proportion[inds] <- proportion
}

# Keep those that are similar with STAT3 in at least half of the cell-lines
dist_filtered <- dist_filtered %>% filter(proportion>=0.5)
GSEAgenes <- unique(c(dist_filtered$cmap_name.x,dist_filtered$cmap_name.y))
GSEAgenes <- GSEAgenes[-which(GSEAgenes=='STAT3')]
GSEAgenes <- data.frame(genes=GSEAgenes)
GSEAgenes <- left_join(GSEAgenes,geneInfo,by=c("genes"="gene_symbol"))
saveRDS(GSEAgenes,'../results/GSEAgenes.rds')