library(tidyverse)

### Load data----

# Load GSEA and latent space candidates
latentgenes <- readRDS('../results/latentgenes.rds')
GSEAgenes <- readRDS('../results/GSEAgenes.rds')

#Load interactions network and paths from source to STAT3
interactions <- read.delim('../preprocessing/preprocessed_data/FilteredOmnipath.tsv') %>% column_to_rownames('X')
paths <- read.csv('../results/pathsToSTAT3.csv') %>% column_to_rownames('X')
paths <- paths %>% filter(paths!='[]')
paths <- paths %>% mutate(pathGenes = str_replace_all(paths, "\\*|\\[|\\]", ""))
paths <- paths %>% mutate(pathGenes = strsplit(pathGenes,",")) %>% unnest(pathGenes)
paths <- paths %>% mutate(pathGenes = gsub("\"", "", pathGenes))
paths <- paths %>% mutate(pathGenes = gsub("\'", "", pathGenes))
paths <- paths %>% mutate(pathGenes=str_replace_all(pathGenes,' ',''))
paths <- paths %>% filter(pathGenes!='')
paths <- paths %>% filter(pathGenes!=' ')
paths <- paths %>% unique()

# First get an interaction subnetwork with all the relevent genes
nodes <- unique(c(paths$source,paths$pathGenes))
interactions <- interactions %>% filter((source %in% nodes) & (target %in% nodes))
nodeAttr <- data.frame(node=unique(c(interactions$source,interactions$target)))
nodeAttr <-  nodeAttr %>% mutate(feature=ifelse(node %in% latentgenes$genes,'latent',
                                                ifelse(node %in% GSEAgenes$genes,'GSEA','other')))
colnames(interactions) <- c("source","sign","target")

# Find now latents in the gene paths to prioritize some candidates from GSEA
paths <- paths %>% filter(pathGenes!='STAT3')
latentGenesInPaths <- left_join(latentgenes,paths,by=c('genes'='pathGenes'))
latentGenesInPaths <- latentGenesInPaths %>% filter(!is.na(paths))
latentGenesInPaths <- latentGenesInPaths %>% select(-source,-paths) %>% unique()

GSEAGenesInPaths <-  left_join(GSEAgenes,paths,by=c('genes'='pathGenes'))
GSEAGenesInPaths <- GSEAGenesInPaths %>% filter(!is.na(paths))
GSEAGenesInPaths <- GSEAGenesInPaths %>% select(-source,-paths) %>% unique()

### Retrive specific geo experiments----
library(limma)
library(GEOquery)
Sys.setenv("VROOM_CONNECTION_SIZE" = 43968200  * 2)
geo <-  getGEO(GEO = 'GSE27869', AnnotGPL = TRUE)
raw <- geo[[1]]

accession <- raw@phenoData@data$geo_accession %>% as.character()
silenced_genes <- raw@phenoData@data$title
silenced_genes <- str_remove_all(silenced_genes,'_KD')

# Keep from original only GSEA candidates and latent candidates for which we have info
GSEAGenesInPaths <- GSEAGenesInPaths %>% filter(genes %in% silenced_genes)
latentGenesInPaths <- latentGenesInPaths %>% filter(genes %in% silenced_genes)

otherGenesInNetwork <- paths %>% filter(!(pathGenes %in% GSEAGenesInPaths$genes)) %>% 
                                 filter(!(pathGenes %in% latentGenesInPaths$genes)) %>%
                                 filter(pathGenes %in% silenced_genes) %>% 
                                 dplyr::select(-source,-paths) %>% unique()


# First perform differential gene expression by utilizing as control STAT3
# The study had only gene siRNA knockdowns not controls.
# Use all the samples at first for statistical significance.

# Get gene expression
exprs_rma <- exprs(raw)

exprs_annotated <- exprs(raw) %>% 
  aggregate(
    by = list(raw@featureData@data$`Gene symbol`), FUN = mean
  ) %>% 
  rename_("Gene symbol" = "Group.1") %>% 
  filter(!grepl("///", `Gene symbol`))
exprs_annotated <- exprs_annotated[2:nrow(exprs_annotated),]
rownames(exprs_annotated) <- NULL
exprs_annotated <- exprs_annotated %>% column_to_rownames("Gene symbol")
exprs_annotated <- as.matrix(exprs_annotated)
exprs_annotated[which(is.na(exprs_annotated))] <- 0
colnames(exprs_annotated) <- silenced_genes
exprs_annotated <- t(scale(t(exprs_annotated)))

# Pathway enrichment analysis for the relevent conditions

releventGenes <- unique(c(GSEAGenesInPaths$genes,latentGenesInPaths$genes,otherGenesInNetwork$pathGenes))
releventGenes <- c(releventGenes,'STAT3')
annotation <- raw@featureData@data %>% dplyr::select('Gene symbol','Gene ID') %>% unique()
annotation <-  annotation %>% filter(`Gene symbol`!='') %>% filter(`Gene symbol`!=' ') %>%
  filter(`Gene ID`!='') %>% filter(`Gene ID`!=' ')
annotation <- annotation %>% filter(`Gene symbol` %in% rownames(exprs_annotated))
# Mannually choose one entrez id for the one with double
annotation <- annotation %>% filter(`Gene ID`!=100009233)
exprs_annotated <- exprs_annotated[annotation$`Gene symbol`,]
print(all(rownames(exprs_annotated)==annotation$`Gene symbol`))
rownames(exprs_annotated) <- annotation$`Gene ID`

# Run fgsea
keegEnrichResults <-  fastenrichment(releventGenes,
                                     rownames(exprs_annotated),
                                     exprs_annotated,
                                     enrichment_space = 'kegg',
                                     n_permutations=10000,
                                     pval_adjustment=F)
keggNES <- keegEnrichResults$NES$`NES KEGG`
keggpadj <- keegEnrichResults$Pval$`Pval KEGG`


# GetSignificants <- function(sampleName,ScoresMatrix,pvalMatrix,p.adj.cutoff=0.01,nes_cutoff=2){
#   pmat <- pvalMatrix[,sampleName]
#   inds <- which(pmat<p.adj.cutoff)
#   smat <- ScoresMatrix[inds,sampleName]
#   inds2 <- which(abs(smat)>nes_cutoff)
#   smat <- smat[inds2]
#   return(smat)
# }
# library(factoextra)
# 
# #Run pca
# pca.samples <- prcomp(t(keggNES),scale=F)
# 
# #Visualize explained variance from top 20 PCs
# png('../figures/scree_geo_validation.png',width=9,height=6,units = "in",res=300)
# fviz_eig(pca.samples,ncp=50)
# dev.off()
# 
# #Built data frame with PCs, labels etc
# df_pca <- pca.samples$x[,1:2]
# df_pca <- data.frame(df_pca)
# df_pca$shRNA <- factor(colnames(keggNES),levels = colnames(keggNES))
# df_pca <- df_pca %>% mutate(space=ifelse(shRNA %in% GSEAgenes$genes,'Kegg GSEA candidate',
#                                          ifelse(shRNA %in% latentgenes$genes,'latent space candidate',
#                                                 'Included in subnetwork')))
# df_pca$space[which(df_pca$shRNA=='STAT3')] <- 'STAT3'
# df_pca$space <- factor(df_pca$space,levels = c('Kegg GSEA candidate',
#                                             'latent space candidate',
#                                             'Included in subnetwork',
#                                             'STAT3'))
# 
# #PCA plot
# pca_plot <- ggplot(df_pca,aes(x=PC1,y=PC2))+ 
#   ggtitle('PCA in KEGG space for GEO validation')
# pca_plot <- pca_plot+geom_point(aes(col =factor(space)),size=3) +
#   scale_color_manual(values = c("Kegg GSEA candidate"="#1ca31c",
#                                 "latent space candidate"="#ff9900",
#                                 "Included in subnetwork"="#BFBFBF",
#                                 'STAT3'='#ff9999')) + 
#   labs(col='')+
#   theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
# print(pca_plot)
#df_kegg <- gather(as.data.frame(keggNES) %>% rownames_to_column('Pathways'),
#                  'Samples','NES',-Pathways)
  
