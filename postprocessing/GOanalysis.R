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

df_annotation = data.frame(gene_id=rownames(cmap))
geneInfo$gene_id <- as.character(geneInfo$gene_id)
df_annotation <- left_join(df_annotation,geneInfo)
rownames(cmap) <- df_annotation$gene_symbol

### Pathway analysis----
# load pathway data
egsea.data(species = "human",returnInfo = TRUE)
print(all(rownames(cmap)==df_annotation$gene_symbol))
rownames(cmap) <- df_annotation$gene_id
keegEnrichResults <-  fastenrichment(sigInfo$sig_id,
                                     geneInfo$gene_id,
                                     cmap,
                                     enrichment_space = 'go_bp',
                                     n_permutations=5000,
                                     pval_adjustment=F)
keggNES <- keegEnrichResults$NES$`NES KEGG`
keggpadj <- keegEnrichResults$Pval$`Pval KEGG`

saveRDS(keggNES,'../results/GONES_STAT3.rds')
saveRDS(keggpadj,'../results/GOpadj_STAT3.rds')