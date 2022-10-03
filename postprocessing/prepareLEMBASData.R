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
sigInfo <- sigInfo %>% mutate(pert_itime= paste(pert_time,pert_time_unit))

# Create a proxy for quality of replicates
# Keep only samples with at least 3 replicates and that satisfy specific conditions.
# Check the LINCS2020 Release Metadata Field Definitions.xlsx file for 
# a complete description of each argument. It can be accessed online 
# or in the data folder.

sigInfo <- sigInfo %>% 
  mutate(quality_replicates = ifelse(is_exemplar_sig==1 & qc_pass==1 & nsample>=3,1,0))
# Filter only experiments of shrnas
sigInfo <- sigInfo %>% filter(pert_type=='trt_sh')

### Plot number of data points per time point
### Boxplot of tas + pvalue per time point
p <- ggplot(sigInfo %>% group_by(pert_itime) %>% mutate(n_points=n()) %>% ungroup() %>%
              mutate(n_points=as.numeric(n_points)) %>%
              mutate(pert_itime=factor(pert_itime,
                                       levels = c('6 h','24 h','48 h','72 h',
                                                             '96 h','120 h','144 h','168 h'))) %>% 
              dplyr::select(pert_itime,n_points) %>% unique(),
            aes(pert_itime,n_points)) +
  geom_bar(stat="identity") + scale_y_log10() +
  xlab('Perturbation time duration (h)') + ylab('Number of data points (log10 scale)')
print(p)

p <- ggplot(sigInfo %>% group_by(pert_itime) %>% mutate(n_points=n()) %>% ungroup() %>%
              mutate(n_points=as.numeric(n_points)) %>%
              mutate(pert_itime=factor(pert_itime,
                                       levels = c('6 h','24 h','48 h','72 h',
                                                  '96 h','120 h','144 h','168 h'))) %>% 
              dplyr::select(pert_itime,tas) %>% unique(),
            aes(pert_itime,tas)) + geom_boxplot()+
  xlab('Perturbation time duration (h)') + ylab('Transcriptional Activity Score')
print(p)


sigInfo <- sigInfo %>% filter(quality_replicates==1)
p <- ggplot(sigInfo %>% group_by(pert_itime) %>% mutate(n_points=n()) %>% ungroup() %>%
              mutate(n_points=as.numeric(n_points)) %>%
              mutate(pert_itime=factor(pert_itime,
                                       levels = c('6 h','24 h','48 h','72 h',
                                                  '96 h','120 h','144 h','168 h'))) %>% 
              dplyr::select(pert_itime,n_points) %>% unique(),
            aes(pert_itime,n_points)) +
  geom_bar(stat="identity") + scale_y_log10() +
  xlab('Perturbation time duration (h)') + ylab('Number of data points (log10 scale)')
print(p)

p <- ggplot(sigInfo %>% group_by(pert_itime) %>% mutate(n_points=n()) %>% ungroup() %>%
              mutate(n_points=as.numeric(n_points)) %>%
              mutate(pert_itime=factor(pert_itime,
                                       levels = c('6 h','24 h','48 h','72 h',
                                                  '96 h','120 h','144 h','168 h'))) %>% 
              dplyr::select(pert_itime,tas) %>% unique(),
            aes(pert_itime,tas)) + geom_boxplot()+
  xlab('Perturbation time duration (h)') + ylab('Transcriptional Activity Score')
print(p)

# Create identifier to signify duplicate
# signatures: meaning same drug, same dose,
# same time duration, same cell-type
sigInfo <- sigInfo %>% 
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))

sigInfo <- sigInfo %>% group_by(duplIdentifier) %>%
  mutate(dupl_counts = n()) %>% ungroup()
gc()
### Find out which shRNAs actually worked
# Load all cmap
#cmap <- readRDS('../../../../../L1000_2021_11_23/cmap_gex_shrna_q1replicates.rds')
#gc()
#cmap <- cmap[,which(colnames(cmap) %in% unique(sigInfo$sig_id))]
#gc()
#cmap_ranked <- apply(cmap,2,rank)
#gc()
#saveRDS(cmap_ranked,'../../../../../L1000_2021_11_23/cmap_gex_shrna_q1replicates_ranked.rds')
#df_annotation = data.frame(gene_id=rownames(cmap))
#geneInfo$gene_id <- as.character(geneInfo$gene_id)
#df_annotation <- left_join(df_annotation,geneInfo)
#print(all(rownames(cmap)==df_annotation$gene_id))
#rownames(cmap) <- df_annotation$gene_symbol

shRNA_kd_significance <- function(sig,sigInfo,cmap){
  gene <- sigInfo$cmap_name[which(sigInfo$sig_id==sig)]
  cell <- sigInfo$cell_iname[which(sigInfo$sig_id==sig)]
  if (gene %in% rownames(cmap)){
    value <- cmap[gene,sig]
    data <- sigInfo %>% filter(cell_iname==cell)
    null_cell_distribution <- cmap[gene,data$sig_id]
    pval <- sum(null_cell_distribution<=value)/length(null_cell_distribution)
  }else{
    pval <- NA
  }
  return(pval)
}

#p.values <- sapply(sigInfo$sig_id,shRNA_kd_significance,sigInfo,cmap)
#print(all(names(p.values)==sigInfo$sig_id))
#saveRDS(p.values,'../results/shRNAs_significance.rds')
p.values <- readRDS('../results/shRNAs_significance.rds')
sigInfo$p.value <- p.values
sigInfo <- sigInfo %>% filter(!is.na(p.value)) %>% filter(p.value<0.05)

### Plot number of data points per time point
### Boxplot of tas + pvalue per time point
p <- ggplot(sigInfo %>% group_by(pert_itime) %>% mutate(n_points=n()) %>% ungroup() %>%
              mutate(n_points=as.numeric(n_points)) %>%
              mutate(pert_itime=factor(pert_itime,
                                       levels = c('6 h','24 h','48 h','72 h',
                                                  '96 h','120 h','144 h','168 h'))) %>% 
              dplyr::select(pert_itime,n_points) %>% unique(),
            aes(pert_itime,n_points)) +
  geom_bar(stat="identity") + scale_y_log10() +
  xlab('Perturbation time duration (h)') + ylab('Number of data points (log10 scale)')
print(p)

### Filter TAS>=0.3 (seems to be needed) to train a lembas model
sigInfo <- sigInfo %>% filter(tas>=0.3)

### Filter multiple time points
sigInfo <- sigInfo %>% mutate(timefree_id = paste0(cmap_name,'_',pert_idose,'_',cell_iname)) %>% group_by(timefree_id) %>%
  mutate(time_points=n_distinct(pert_itime)) %>% mutate(tas_max=max(tas)) %>%
  mutate(timekeep=ifelse(tas==tas_max | time_points==1,pert_itime,NA)) %>%  
  mutate(timekeep = unique(timekeep)[which(!is.na(unique(timekeep)))]) %>% ungroup() %>%
  filter(pert_itime==timekeep) %>% dplyr::select(-timekeep,-time_points,-timefree_id,-tas_max) %>% unique()


# First load Omnipath prior knowledge network (PKN)
#library(OmnipathR)
#interactions <- import_omnipath_interactions()
#interactions <- interactions %>% mutate(kegg=grepl(pattern="KEGG",x=sources)) %>% 
#  mutate(signoir=grepl(pattern="SIGNOR",x=sources)) %>% filter(kegg==T | signoir==T)
#interactions <- interactions  %>% filter(n_resources>1) %>% dplyr::select(c('source'='source_genesymbol'),
#                                                                          c('target'='target_genesymbol'),
#                                                                          is_inhibition,is_stimulation) %>% unique()

#interactions <- interactions %>% mutate(interaction=ifelse(is_stimulation==1,
#                                                           ifelse(is_inhibition!=1,1,0),
#                                                           ifelse(is_inhibition!=0,-1,0))) %>%
#  dplyr::select(source,interaction,target) %>% unique()

#write.table(interactions, file = '../preprocessing/preprocessed_data/FilteredOmnipath_v2.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
## Load omnipath
#interactions <- read.delim('../preprocessing/preprocessed_data/FilteredOmnipath_v2.tsv') %>% column_to_rownames('X')
## Filter shRNAs not in OmniPath
#sigInfo <- sigInfo %>% filter(cmap_name %in% unique(c(interactions$source,interactions$target))) %>% unique()

### Load CCLE data----
ccle <- t(data.table::fread('../data/CCLE/CCLE_expression_old.csv') %>% column_to_rownames('V1'))
ccle <- as.data.frame(ccle) %>% rownames_to_column('V1') %>% separate(V1,c('gene_id','useless'),sep=" ") %>%
  dplyr::select(-useless) %>% column_to_rownames('gene_id')
ccle <- as.data.frame(t(ccle)) %>% rownames_to_column('DepMap_ID')
sample_info <- data.table::fread('../data/CCLE/sample_info_old.csv') %>% dplyr::select(DepMap_ID,stripped_cell_line_name) %>%
  unique()
ccle <- left_join(ccle,sample_info) %>% dplyr::select(-DepMap_ID) %>%
  column_to_rownames('stripped_cell_line_name')
ccle <- ccle[which(rownames(ccle) %in% unique(sigInfo$cell_iname)),]
write.table(ccle, file = '../data/CCLE/preprocessed_ccle_tasfiltered.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

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
#TF_activities <- TF_activities[,which(colnames(TF_activities) %in% unique(c(interactions$source,interactions$target)))]
write.table(TF_activities, file = '../results/filtered_shrnas_tf_activities.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)


### Create signaling model------------------
annot <- read.delim('../data/annotation/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab')
annot <- annot %>% dplyr::select(Gene.names...primary..,Entry) %>% unique()
colnames(annot) <- c('cmap_name','cmap_name_uniprot')
annot <- annot %>% filter(cmap_name!='') %>% mutate(cmap_name=strsplit(cmap_name,';')) %>% unnest(cmap_name) %>% unique()

### Targeted TFS
tfs_tageted <- colnames(TF_activities)[which(colnames(TF_activities) %in% unique(sigInfo$cmap_name))]
tfs_tageted <- as.data.frame(tfs_tageted)
colnames(tfs_tageted) <- 'Entry'
write.table(tfs_tageted, file = '../preprocessing/preprocessed_data/targetd_tfs.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

### FILTER SHRNAS AND TFS NOT IN THE TRIMMED NETWORK
#pkn <- read.delim('../preprocessing/preprocessed_data/KEGG-Model.tsv')
#pkn <- read.delim('../preprocessing/preprocessed_data/KEGG_Ligands-Model.tsv')
pkn <- read.delim('../preprocessing/preprocessed_data/macrophageL1000_Ligands-Model.tsv')
omnipath <- read.delim('../data/annotation/omnipath_webservice_interactions__recent.tsv')
human <- 9606
omnipath <- omnipath %>% filter(ncbi_tax_id_source==human & ncbi_tax_id_target==human)
omnipathAnnot <- rbind(omnipath %>% dplyr::select(c("cmap_name"="source_genesymbol"),c("cmap_name_uniprot"="source")),
                       omnipath %>% dplyr::select(c("cmap_name"="target_genesymbol"),c("cmap_name_uniprot"="target"))) %>% 
  unique()
print(length(which(sigInfo$cmap_name %in% omnipathAnnot$cmap_name)))
sigInfo <- left_join(sigInfo,omnipathAnnot)
sigInfo <- sigInfo %>% filter((cmap_name_uniprot %in% unique(c(pkn$source,pkn$target))))%>% 
  filter(!is.na(cmap_name_uniprot)) %>% unique()
tfs <- data.frame("cmap_name"=colnames(TF_activities))
tfs <- left_join(tfs,annot)
print(all(colnames(TF_activities)==tfs$cmap_name))
colnames(TF_activities) <- tfs$cmap_name_uniprot
TF_activities <- TF_activities[,which(colnames(TF_activities) %in% unique(c(pkn$source,pkn$target)))]
write.table(TF_activities, file = '../results/trimmed_shrnas_tf_activities_ligandpkn_old.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

### Filter multiple time points
sigInfo <- sigInfo %>% mutate(timefree_id = paste0(cmap_name,'_',pert_idose,'_',cell_iname)) %>% group_by(timefree_id) %>%
  mutate(time_points=n_distinct(pert_itime)) %>% mutate(tas_max=max(tas)) %>%
  mutate(timekeep=ifelse(tas==tas_max | time_points==1,pert_itime,NA)) %>%  
  mutate(timekeep = unique(timekeep)[which(!is.na(unique(timekeep)))]) %>% ungroup() %>%
  filter(pert_itime==timekeep) %>% dplyr::select(-timekeep,-time_points,-timefree_id,-tas_max) %>% unique()

### FILTER CCLE GENES NOT IN THE NETWORK ---> those not in the trimmed omnipath.
ccleAnnotation <- omnipathAnnot %>% filter(cmap_name %in% colnames(ccle))
ccleAnnotation <- left_join(rbind(pkn %>% dplyr::select(c("cmap_name_uniprot"="source")) %>% unique(),
                                  pkn %>% dplyr::select(c("cmap_name_uniprot"="target")) %>% unique()) %>% unique(),
                            omnipathAnnot)
ccle <- ccle[,which(colnames(ccle) %in% ccleAnnotation$cmap_name)]
ccleAnnotation <-  ccleAnnotation %>% filter(cmap_name %in% colnames(ccle))
#ccleAnnotation %>% group_by(cmap_name) %>% summarise(counts=n_distinct(cmap_name_uniprot)) %>% filter(counts>1)
# Manually choose the name for the name for the doubles
ccleAnnotation$cmap_name_uniprot[which(ccleAnnotation$cmap_name=='CDKN2A')] <- 'P42771'
ccleAnnotation$cmap_name_uniprot[which(ccleAnnotation$cmap_name=='GNAS')] <- 'O95467'
ccleAnnotation <- ccleAnnotation %>% unique()
ccle <- ccle[,ccleAnnotation$cmap_name]
which(is.na(colnames(ccle)))
which(is.na(ccle))
print(all(ccleAnnotation$cmap_name==colnames(ccle)))
colnames(ccle) <- ccleAnnotation$cmap_name_uniprot
ccle <- ccle[,which(colnames(ccle) %in% unique(c(pkn$source,pkn$target)))]
ccle <- ccle[which(rownames(ccle) %in% unique(sigInfo$cell_iname)),]
sigInfo <- sigInfo %>% filter(cell_iname %in% rownames(ccle)) %>% unique()
ccle <- ccle[which(rownames(ccle) %in% unique(sigInfo$cell_iname)),]
#ccle <- ccle[,which(apply(ccle,2,sum)!=0)]
write.table(ccle, file = '../data/CCLE/trimmed_ccle_ligandpkn_old.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

### Ligands modeling--------------------------------------------------------------
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

# keep samples with at least 2 replicates instead of 3
sigInfo <- sigInfo %>% 
  mutate(quality_replicates = ifelse(is_exemplar_sig==1 & qc_pass==1 & nsample>=3,1,0))
# Filter only experiments of ligands
sigInfo <- sigInfo %>% filter(pert_type=='trt_lig')
sigInfo <- sigInfo %>% filter(quality_replicates==1)
sigInfo <- sigInfo %>% filter(tas>0.25)

# Create identifier to signify duplicate
# signatures: meaning same drug, same dose,
# same time duration, same cell-type
sigInfo <- sigInfo %>% 
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))

sigInfo <- sigInfo %>% group_by(duplIdentifier) %>%
  mutate(dupl_counts = n()) %>% ungroup()
gc()

#pkn <- read.delim('../preprocessing/preprocessed_data/KEGG_Ligands-Model.tsv')
pkn <- read.delim('../preprocessing/preprocessed_data/macrophageL1000_Ligands-Model.tsv')
omnipath <- read.delim('../data/annotation/omnipath_webservice_interactions__recent.tsv')
human <- 9606
omnipath <- omnipath %>% filter(ncbi_tax_id_source==human & ncbi_tax_id_target==human)
omnipathAnnot <- rbind(omnipath %>% dplyr::select(c("cmap_name"="source_genesymbol"),c("cmap_name_uniprot"="source")),
                       omnipath %>% dplyr::select(c("cmap_name"="target_genesymbol"),c("cmap_name_uniprot"="target"))) %>% 
  unique()
print(length(which(sigInfo$cmap_name %in% omnipathAnnot$cmap_name)))
sigInfo <- left_join(sigInfo,omnipathAnnot)
sigInfo <- sigInfo %>% filter((cmap_name_uniprot %in% pkn$source))

### Filter multiple time points
sigInfo <- sigInfo %>% filter(pert_time<24)
sigInfo <- sigInfo %>% mutate(timefree_id = paste0(cmap_name,'_',pert_idose,'_',cell_iname)) %>% group_by(timefree_id) %>%
  mutate(time_points=n_distinct(pert_itime)) %>% mutate(tas_max=max(tas)) %>%
  mutate(timekeep=ifelse(tas==tas_max | time_points==1,pert_itime,NA)) %>%  
  mutate(timekeep = unique(timekeep)[which(!is.na(unique(timekeep)))]) %>% ungroup() %>%
  filter(pert_itime==timekeep) %>% dplyr::select(-timekeep,-time_points,-timefree_id,-tas_max) %>% unique()

# Filter cell-lines not in ccle
### Load CCLE data
ccle <- t(data.table::fread('../data/CCLE/CCLE_expression_old.csv') %>% column_to_rownames('V1'))
ccle <- as.data.frame(ccle) %>% rownames_to_column('V1') %>% separate(V1,c('gene_id','useless'),sep=" ") %>%
  dplyr::select(-useless) %>% column_to_rownames('gene_id')
ccle <- as.data.frame(t(ccle)) %>% rownames_to_column('DepMap_ID')
sample_info <- data.table::fread('../data/CCLE/sample_info_old.csv') %>% dplyr::select(DepMap_ID,stripped_cell_line_name) %>%
  unique()
ccle <- left_join(ccle,sample_info) %>% dplyr::select(-DepMap_ID) %>%
  column_to_rownames('stripped_cell_line_name')

ccleAnnotation <- omnipathAnnot %>% filter(cmap_name %in% colnames(ccle))
ccleAnnotation <- left_join(rbind(pkn %>% dplyr::select(c("cmap_name_uniprot"="source")) %>% unique(),
                                  pkn %>% dplyr::select(c("cmap_name_uniprot"="target")) %>% unique()) %>% unique(),
                            omnipathAnnot)
ccle <- ccle[,which(colnames(ccle) %in% ccleAnnotation$cmap_name)]
ccleAnnotation <-  ccleAnnotation %>% filter(cmap_name %in% colnames(ccle))
#ccleAnnotation %>% group_by(cmap_name) %>% summarise(counts=n_distinct(cmap_name_uniprot)) %>% filter(counts>1)
# Manually choose the name for the name for the doubles
ccleAnnotation$cmap_name_uniprot[which(ccleAnnotation$cmap_name=='CDKN2A')] <- 'P42771'
ccleAnnotation$cmap_name_uniprot[which(ccleAnnotation$cmap_name=='GNAS')] <- 'O95467'
ccleAnnotation <- ccleAnnotation %>% unique()
ccle <- ccle[,ccleAnnotation$cmap_name]
which(is.na(colnames(ccle)))
which(is.na(ccle))
print(all(ccleAnnotation$cmap_name==colnames(ccle)))
colnames(ccle) <- ccleAnnotation$cmap_name_uniprot
ccle <- ccle[,which(colnames(ccle) %in% unique(c(pkn$source,pkn$target)))]
ccle <- ccle[which(rownames(ccle) %in% unique(sigInfo$cell_iname)),]
sigInfo <- sigInfo %>% filter(cell_iname %in% rownames(ccle)) %>% unique()
ccle <- ccle[which(rownames(ccle) %in% unique(sigInfo$cell_iname)),]
#ccle <- ccle[,which(apply(ccle,2,sum)!=0)]
write.table(ccle, file = '../data/CCLE/trimmed_ccle_ligands_old.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

### Get rid of multiple duplicates (keep only one for now)
sigInfo <- sigInfo %>% 
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname)) %>% group_by(duplIdentifier) %>%
  mutate(dupl_counts = n()) %>% ungroup()
sigInfo$dupl_keep <- FALSE
for (i in 1:nrow(sigInfo)){
  if (sigInfo$dupl_counts[i]>1){
    tmp <- sigInfo %>% filter(duplIdentifier==sigInfo$duplIdentifier[i])
    max_tas = max(tmp$tas)
    sig <- tmp$sig_id[which(tmp$tas==max_tas)]
    sigInfo$dupl_keep[which(sigInfo$sig_id==sig)] <- TRUE
  }else{
    sigInfo$dupl_keep[i] <- TRUE
  }
  if (i%%100==0){
    print(paste0('Finished:',i))
  }
}
sigInfo <- sigInfo %>% filter(dupl_keep==TRUE) %>% dplyr::select(-dupl_keep) %>% unique()

### Infer transcription factors with Dorothea
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
ds_path <- '../../../../../L1000_2021_11_23/level5_beta_trt_misc_n8283x12328.gctx'
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
annot <- read.delim('../data/annotation/uniprot-reviewed_yes+AND+organism__Homo+sapiens+(Human)+[9606]_.tab')
annot <- annot %>% dplyr::select(Gene.names...primary..,Entry) %>% unique()
colnames(annot) <- c('cmap_name','cmap_name_uniprot')
tfs <- data.frame("cmap_name"=colnames(TF_activities))
tfs <- left_join(tfs,annot)
print(all(colnames(TF_activities)==tfs$cmap_name))
colnames(TF_activities) <- tfs$cmap_name_uniprot
TF_activities <- TF_activities[,which(colnames(TF_activities) %in% unique(c(pkn$source,pkn$target)))]
write.table(TF_activities, file = '../results/trimmed_ligands_old_tf_activities.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)


### Get condition matrix (log dose)
print(unique(sigInfo$pert_dose_unit))
min_val <- 0.01 # ng/mL
sigInfo <- sigInfo %>% mutate(log10Dose = ifelse(pert_dose_unit=="ng/uL",log10(1000*pert_dose/min_val+1),log10(pert_dose/min_val+1)))
hist(sigInfo$log10Dose)

conditionMatrix <- sigInfo %>% dplyr::select(sig_id,cmap_name_uniprot,log10Dose) %>% unique()
conditionMatrix <- as.matrix(conditionMatrix %>% spread(cmap_name_uniprot,log10Dose) %>% column_to_rownames('sig_id'))
conditionMatrix[which(is.na(conditionMatrix))] <- 0.0
write.table(conditionMatrix,'../results/Ligands_conditions_old.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

cellInfo <- sigInfo %>% dplyr::select(sig_id,cell_iname) %>% unique()
cellInfo <- cellInfo %>% mutate(value=1) %>% spread('cell_iname','value')
cellInfo[is.na(cellInfo)] <- 0
write.table(cellInfo, file = '../preprocessing/preprocessed_data/all_filtered_cells_ligand_old.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

ligands <- colnames(conditionMatrix)

### Split random 10fold validation----

### Only for cell-line specific models
#sigInfo <- sigInfo %>% group_by(cell_iname) %>% mutate(samples_per_cell=n_distinct(sig_id)) %>% ungroup()
#sigInfo <- sigInfo %>% filter(cell_iname=='HT29')

### Get rid of multiple duplicates (keep only one for now)
sigInfo <- sigInfo %>% 
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname)) %>% group_by(duplIdentifier) %>%
  mutate(dupl_counts = n()) %>% ungroup()
sigInfo$dupl_keep <- FALSE
for (i in 1:nrow(sigInfo)){
  if (sigInfo$dupl_counts[i]>1){
    tmp <- sigInfo %>% filter(duplIdentifier==sigInfo$duplIdentifier[i])
    max_tas = max(tmp$tas)
    sig <- tmp$sig_id[which(tmp$tas==max_tas)]
    sigInfo$dupl_keep[which(sigInfo$sig_id==sig)] <- TRUE
  }else{
    sigInfo$dupl_keep[i] <- TRUE
  }
  if (i%%100==0){
    print(paste0('Finished:',i))
  }
}
sigInfo <- sigInfo %>% filter(dupl_keep==TRUE) %>% dplyr::select(-dupl_keep) %>% unique()

library(caret)
total_samples <- sigInfo$sig_id
folds <- createFolds(total_samples, k = 10, list = TRUE, returnTrain = TRUE)

i <- 0
for (fold in folds){
  samples <- total_samples[fold]
  
  train_samples <- sigInfo %>% filter(sig_id %in% samples)
  val_samples <- sigInfo %>% filter(!(sig_id %in% samples))
  
  data.table::fwrite(as.data.frame(train_samples),paste0('../data/10fold_cross_validation/LigandPKN_old/train_sample_',i,'.csv'),row.names = T)
  data.table::fwrite(as.data.frame(val_samples),paste0('../data/10fold_cross_validation/LigandPKN_old/val_sample_',i,'.csv'),row.names = T)
  
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
  
  data.table::fwrite(as.data.frame(train_samples),paste0('../data/3fold_cross_validation/cell_based/LigandPKN_old/train_sample_',i,'.csv'),row.names = T)
  data.table::fwrite(as.data.frame(val_samples),paste0('../data/3fold_cross_validation/cell_based/LigandPKN_old/val_sample_',i,'.csv'),row.names = T)
  
  i <- i+1
}

### Knock-outs are simulated as putting -5 in the input node----
# Create a general condition input data (there are not multiple doses (only NAs and a specific dose))
# Probably finding the cmap_name while training and saying -5 there is sufficient
sigInfo <- sigInfo %>% filter(cmap_name_uniprot %in% unique(c(pkn$source,pkn$target)))
cellInfo <- sigInfo %>% dplyr::select(sig_id,cell_iname) %>% unique()
sigInfo <- sigInfo %>% dplyr::select(sig_id,cmap_name_uniprot) %>% unique()
sigInfo <- sigInfo %>% mutate(value=-10) %>% spread('cmap_name_uniprot','value')
sigInfo[is.na(sigInfo)] <- 0
write.table(sigInfo, file = '../preprocessing/preprocessed_data/all_filtered_Kds_ligandpkn_old.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
cellInfo <- cellInfo %>% mutate(value=1) %>% spread('cell_iname','value')
cellInfo[is.na(cellInfo)] <- 0
cellInfo <- cellInfo %>% column_to_rownames('sig_id')
print(all(colnames(cellInfo)==rownames(ccle)))
cellInfo <- cellInfo[,rownames(ccle)]
cellInfo <- cellInfo %>% rownames_to_column('sig_id')
write.table(cellInfo, file = '../preprocessing/preprocessed_data/all_filtered_cells_ligandpkn_old.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
