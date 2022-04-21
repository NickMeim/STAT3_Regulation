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

### Keep only STAT3 knockdowns
sigInfo <- sigInfo %>% filter(cmap_name=='STAT3')

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

col_fun = colorRamp2(c(-4, 0, 4), c("red", "white", "blue"))
png('../figures/STAT3_GeX.png',width=16,height=8,units = "in",res=300)
Heatmap(cmap, col=col_fun,
        row_title = 'Genes',
        column_title = 'Samples',
        show_row_names =F,
        show_column_names  =F,
        heatmap_legend_param = list(title = 'z-score'))
dev.off()
df_annotation = data.frame(gene_id=rownames(cmap))
geneInfo$gene_id <- as.character(geneInfo$gene_id)
df_annotation <- left_join(df_annotation,geneInfo)
rownames(cmap) <- df_annotation$gene_symbol

### Infer transcription factors with Dorothea----

minNrOfGenes = 5

# Load requeired packages
library("dorothea")
dorotheaData = read.table('../data/dorothea.tsv', sep = "\t", header=TRUE)
confidenceFilter = is.element(dorotheaData$confidence, c('A', 'B'))
dorotheaData = dorotheaData[confidenceFilter,]

# Estimate TF activities
settings = list(verbose = TRUE, minsize = minNrOfGenes)
TF_activities = run_viper(cmap, dorotheaData, options =  settings)

#TF_activities <- t(TF_activities)
write.table(TF_activities, file = '../results/stat3_shrna_tf_activities.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)
#tf_ranks <- apply(TF_activities,2,rank)
#tf_median_rank <- as.matrix(apply(tf_ranks,1,median))
#tf_median_rank[which(rownames(tf_median_rank)=='STAT3')]

### Run CARNIVAL to infer signalling network----
library(CARNIVAL)
# First load Omnipath prior knowledge network (PKN)
library(OmnipathR)
interactions <- import_omnipath_interactions()
interactions <- interactions  %>% filter(n_resources>1) %>% dplyr::select(c('source'='source_genesymbol'),
                                                                   c('target'='target_genesymbol'),
                                                                   is_inhibition,is_stimulation) %>% unique()

interactions <- interactions %>% filter(!(is_inhibition==0 & is_stimulation==0)) %>% unique()
interactions <- interactions %>% mutate(interaction=ifelse(is_stimulation!=0,1,-1)) %>%
  dplyr::select(source,interaction,target) %>% unique()
#write.table(interactions, file = '../preprocessing/preprocessed_data/FilteredOmnipath.tsv', quote=FALSE, sep = "\t", row.names = TRUE, col.names = NA)

# Get top-bottom tfs
top_bot_indices <- function(v,num){
  top <- order(v,decreasing = TRUE)[1:num]
  bot <-  order(v)[1:num]
  return(c(top,bot))
}


log_con <- file("log.txt", open="a")
for (j in 1:ncol(TF_activities)){
  
  cat(paste0('Iteration ',j,'/',ncol(TF_activities)), file = log_con, sep="\n")
  tf_activities <- TF_activities[which(rownames(TF_activities) %in% c(interactions$source,
                                                                      interactions$target)),]
  ind <- top_bot_indices(tf_activities[,j],30)
  tf_activities <- tf_activities[ind,]
  tf_activities <- tf_activities %>% dplyr::select(colnames(tf_activities)[j])
  names <- rownames(tf_activities)
  tf_activities <- tf_activities[,1]
  names(tf_activities) <- names
  
  # Run carnival
  # YOU MUST FIND WHERE CPLEX IS INSTALLED IN YOUR OWN COMPUTER
  CplexPath <- 'C:/Program Files/IBM/ILOG/CPLEX_Studio201/cplex/bin/x64_win64/cplex.exe'
  carnivalOptions <- defaultCplexCarnivalOptions()
  carnivalOptions$solverPath <- CplexPath
  carnivalOptions$timelimit <- 1200
  # Output dir
  Result_dir <- paste0("../results/stat3_networks/",colnames(TF_activities)[j])
  dir.create(Result_dir, showWarnings = FALSE)
  carnivalOptions$outputFolder <- Result_dir
  
  inverseCarnivalResults <- runInverseCarnival( measurements = tf_activities, 
                                                priorKnowledgeNetwork = interactions, 
                                                carnivalOptions = carnivalOptions)
  
  # Save interaction networks
  nets <- inverseCarnivalResults$sifAll
  nodes <- inverseCarnivalResults$attributesAll
  for (i in 1:length(nets)){
    t <- nets[[i]]
    t <- as.data.frame(t)
    t$Node1 <- as.character(t$Node1)
    t$Node2 <- as.character(t$Node2)
    t$Sign <- as.character(t$Sign)
    write_tsv(t,paste0(Result_dir,'/','interactions_1_model',i,'.tsv'))
    t <- nodes[[i]]
    t <- as.data.frame(t)
    t$Nodes <- as.character(t$Nodes)
    t$Activity <- as.numeric(t$Activity)
    write_delim(t,paste0(Result_dir,'/','nodesActivity_1_model',i,'.txt'),delim = '\t')
  }
  t <- as.data.frame(inverseCarnivalResults$weightedSIF)
  t$Node1 <- as.character(t$Node1)
  t$Node2 <- as.character(t$Node2)
  t$Sign <- as.character(t$Sign)
  t$Weight <- as.numeric(t$Weight)
  write_delim(t,paste0(Result_dir,'/','weightedModel_1.txt'),delim = '\t')
  t <- as.data.frame(inverseCarnivalResults[["nodesAttributes"]])
  write_delim(t,paste0(Result_dir,'/','nodesAttributes_1.txt'),delim = '\t')
}
close(log_con)


## Infer a weighted representative network for STAT3
net_files <- list.files(path='../results/stat3_networks',recursive = T,full.names = T) 
net_files <- as.data.frame(net_files)
net_files <- net_files %>%
  mutate(meas = grepl(pattern = "meas_",x = net_files),
         log = grepl(pattern = ".log",x = net_files),
         time = grepl(pattern = "elapsed_time.txt",x = net_files),
         res = grepl(pattern = "results_CARNIVAL.Rdata",x = net_files),
         empty = grepl(pattern = "emptyNetwork",x = net_files),
         lp = grepl(pattern = 'lpFile',x=net_files))
net_files <-  net_files %>% filter(meas == F & log == F & time == F & res == F & empty == F & lp==F)

match_weight <- c("nodesAttributes_1.txt","weightedModel_1.txt")
net_files <- net_files %>% dplyr::select(net_files) %>% 
  mutate(weighted = grepl(pattern = paste(match_weight,collapse = "|"),x = net_files))
net_files <- net_files %>% filter(weighted == F)

# Check which are edges and which are node attributes
net_files <- net_files %>% mutate(edgeInfo = grepl(pattern = "interactions",x = net_files),
                                  nodeInfo = grepl(pattern = "nodesActivity",x = net_files)) %>%
  dplyr::select(-weighted)
nodeFiles <- net_files %>% filter(nodeInfo==T)
edgeFiles <- net_files %>% filter(edgeInfo==T)

nodes <- data.frame()
edges <- data.frame()
for (i in 1:nrow(nodeFiles)){
  nodes <- rbind(nodes,
                 read.delim(nodeFiles$net_files[i]))
  edges <- rbind(edges,
                 read.delim(edgeFiles$net_files[i]))
  message('Processing image ', i, ' of ', nrow(nodeFiles))
}

nodes <- aggregate(Activity~Nodes,data=nodes,mean)
nodes <- nodes %>% mutate(AbsActivity=abs(Activity))
nodes <- nodes %>% filter(AbsActivity>0.5)
edges <- aggregate(Sign~Node1+Node2,data=edges,mean)
edges <- edges %>% mutate(weight=abs(Sign),
                          Sign=sign(Sign))
edges <- edges %>% filter(weight>0.5) %>% filter((Node1 %in% nodes$Nodes) & (Node2 %in% nodes$Nodes))
nodes <- nodes %>% filter(Nodes %in% unique(c(edges$Node1,edges$Node1)))

### Pathway analysis----
# load pathway data
egsea.data(species = "human",returnInfo = TRUE)
print(all(rownames(cmap)==df_annotation$gene_symbol))
rownames(cmap) <- df_annotation$gene_id
keegEnrichResults <-  fastenrichment(sigInfo$sig_id,
                                     geneInfo$gene_id,
                                     cmap,
                                     enrichment_space = 'kegg',
                                     n_permutations=5000)
keggNES <- keegEnrichResults$NES$`NES KEGG`
keggpadj <- keegEnrichResults$Pval$`Pval KEGG`

saveRDS(keggNES,'../results/keggNES_STAT3.rds')
saveRDS(keggpadj,'../results/keggpadj_STAT3.rds')

# Get significantly regulated pathways
# get_signigicant <- function(pvals,level=0.05){
#   return(which(pvals<level))
# }
get_signigicant_counts <- function(pvals,level=0.05){
  return(length(which(pvals<level))/length(pvals))
}
significants <- apply(keggpadj, 1 ,get_signigicant_counts,level=0.05)
significants <- which(significants>=0.5)
#significants <- Reduce(intersect,significants)

hist(keggNES[significants,],breaks=20)
keggNES_signigicants <- keggNES[significants,]

df_kegg <- gather(as.data.frame(keggNES_signigicants) %>% rownames_to_column('Pathways'),
                  'Samples','NES',-Pathways)
png(file="../figures/keggHeatmap.png",width=16,height=8,units = "in",res=600)
ggplot(df_kegg, aes(Samples,Pathways, fill= NES)) + 
  geom_tile()+theme(axis.text.x=element_blank(),text = element_text(size=13))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_gradient2()+ ggtitle("STAT3 regulated pathways")

dev.off()

keggNES_signigicants_mean <- as.matrix(apply(keggNES_signigicants,1,mean))

df_kegg_mean <- as.data.frame(keggNES_signigicants_mean) %>% rownames_to_column('Pathways')
png(file="../figures/MeankeggBarplot.png",width=8,height=8,units = "in",res=600)
ggplot(df_kegg_mean, aes(V1,Pathways, fill= V1)) + ylab('')+xlab('')+
  guides(fill=guide_legend(title="NES"))+
  geom_bar(stat='identity')+
  theme_classic()+
  theme(text = element_text(size=13),
        plot.title = element_text(hjust = 0.5),
        legend.position="none")+
  scale_fill_gradient2()+ ggtitle("Average NES")

dev.off()

# Find the downregulated ones and get their genes as potential targets
# Compare their available shRNAs with STAT3 shRNA
pathways_down <- rownames(keggNES_signigicants_mean)[which(keggNES_signigicants_mean<(-0.3))]

split_geneset_id <- function(name,pattern){
  new_name <- str_split(name,pattern)[[1]][2]
  return(new_name)
}
pathways_down <- sapply(pathways_down,split_geneset_id , pattern = 'FL1000_KEGG_')

human_paths <- kegg.pathways$human$kg.sets
gene_candidates <- human_paths[pathways_down]
gene_candidates <- Reduce(union,gene_candidates)

pathways_up <- rownames(keggNES_signigicants_mean)[which(keggNES_signigicants_mean>0)]
pathways_up <- sapply(pathways_up,split_geneset_id , pattern = 'FL1000_KEGG_')
human_paths <- kegg.pathways$human$kg.sets
gene_exclude  <- human_paths[pathways_up]
gene_exclude <- Reduce(union,gene_exclude)

gene_candidates <- gene_candidates[which(!(gene_candidates %in% gene_exclude))]
gene_candidates <- gene_candidates[which(gene_candidates %in% geneInfo$gene_id)]
df_gene_canditates <- data.frame(gene_id=as.numeric(gene_candidates))
df_gene_canditates <- left_join(df_gene_canditates,geneInfo)
saveRDS(df_gene_canditates,'../results/gene_candidates.rds')

### Therapeutic screening----

# Screen drugs and shRNAs based on the GSEA distance with the STAT3 knockdown
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

# Load gene candidates
df_gene_canditates <- readRDS('../results/gene_candidates.rds')

# Keep L1000 shRNA data of candidate genes
sigInfo <- sigInfo %>% filter(cmap_name %in% df_gene_canditates$gene_symbol)

# Now exclude low TAS (signal) signatures cause you are going to compare with STAT3
sigInfo <- sigInfo %>% filter(tas>=0.3)

# Calculate distance of candidates from STAT3
keggNES <- readRDS('../results/keggNES_shRNAs.rds')
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

## Find neighbors from autoencoder----
# Load embeddings
emb <- read.csv('../../shrna_embs1024.csv')
#emb <- emb %>% mutate(shRNA='other')
#emb$shRNA[which(emb$X %in% stat3$sig_id)] <- 'STAT3'
emb <- emb %>% column_to_rownames('X')
emb <- as.matrix(emb)
gc()

sigInfo <- read.delim('../data/siginfo_beta.txt')
sigInfo <- sigInfo %>% 
  mutate(quality_replicates = ifelse(is_exemplar_sig==1 & qc_pass==1 & nsample>=3,1,0))
sigInfo <- sigInfo %>% filter(pert_type=='trt_sh')
sigInfo <- sigInfo %>% filter(quality_replicates==1)
sigInfo <- sigInfo %>% 
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))

sigInfo <- sigInfo %>% group_by(duplIdentifier) %>%
  mutate(dupl_counts = n()) %>% ungroup()
stat3Info <- sigInfo %>% filter(cmap_name=='STAT3')
#sigInfo <- sigInfo %>% filter(tas>=0.3)
#sigInfo <- rbind(sigInfo,stat3Info)
#sigInfo <- sigInfo %>% unique()

emb <- emb[which(rownames(emb) %in% sigInfo$sig_id),] 
#emb <- emb[sigInfo$sig_id,]

# # Euclidean distances or cosine similarities in latent space
# #distance <- as.matrix(dist(emb, method = "euclidean",diag = F,upper = F))
# library(lsa)
# X <- t(emb)
# distance <- cosine(X)
# colnames(distance) <- rownames(emb)
# rownames(distance) <- rownames(emb)
# 
# ### Convert matrix into data frame
# # Keep only unique (non-self) pairs
# distance[lower.tri(distance,diag = T)] <- -100
# dist <- reshape2::melt(distance)
# dist <- dist %>% filter(value != -100)
# 
# sigInfo <- sigInfo %>% select(sig_id,cmap_name,cell_iname,pert_idose,pert_itime,duplIdentifier) %>% unique()
# 
# # Merge meta-data info and distances values
# dist <- left_join(dist,sigInfo,by = c("Var1"="sig_id"))
# dist <- left_join(dist,sigInfo,by = c("Var2"="sig_id"))
# dist <- dist %>% mutate(is_duplicate = (duplIdentifier.x==duplIdentifier.y))
# dist <- dist %>% filter(!is.na(value))
# #saveRDS(dist,'../results/sim_cosine_latent_5dupls.rds')
# gc()
# 
# # Keep only STAT3 with other candidate pairs
# dist <- dist %>% filter(cmap_name.x=='STAT3' | cmap_name.y=='STAT3')
# dist <- dist %>% filter(!(cmap_name.x=='STAT3' & cmap_name.y=='STAT3'))
# dist <- dist %>% filter(cell_iname.x==cell_iname.y)
# dist <- dist %>% mutate(cell_iname = cell_iname.x) %>% dplyr::select(-cell_iname.x,-cell_iname.y) %>% unique()
# 
# # # First aggregate duplicates
# # dist <- dist %>% mutate(pairID=ifelse(cmap_name.x!='STAT3',duplIdentifier.x,duplIdentifier.y)) %>%
# #   group_by(pairID) %>% mutate(med_value=median(value)) %>% ungroup()
# # dist$value <- dist$med_value
# # dist <- dist %>% dplyr::select(-med_value,-Var1,-Var2,-duplIdentifier.x,-duplIdentifier.y) %>% unique()
# #df <- dist %>% group_by(pairID) %>% summarise(n())
# 
# # Get neighbors per cell
# cells <- unique(as.character(dist$cell_iname))
# 
# # Filter dist
# dist_filtered <- dist %>% filter(value>=0.3) #from duplicates distr and because it is cosine
# 
# # Get counts of how many times each candidate was found as a neighbor
# genes <- unique(c(dist_filtered$cmap_name.x,dist_filtered$cmap_name.y))
# genes <- genes[-which(genes=='STAT3')]
# 
# dist_filtered <- dist_filtered %>% mutate(proportion=0)
# for (i in 1:length(genes)){
#   inds <- which(dist_filtered$cmap_name.x==genes[i] | dist_filtered$cmap_name.y==genes[i])
#   df <- dist %>% filter(cmap_name.x==genes[i] | cmap_name.y==genes[i]) %>% dplyr::select(-cmap_name.x,-cmap_name.y) %>% unique()
#   df_filtered <- dist_filtered %>% filter(cmap_name.x==genes[i] | cmap_name.y==genes[i]) %>% dplyr::select(-cmap_name.x,-cmap_name.y) %>% unique()
#   proportion <- length(unique(df_filtered$cell_iname))/length(unique(df$cell_iname))
#   dist_filtered$proportion[inds] <- proportion
# }

### After clustering anlysis get 250 clusters
km <- kmeans(emb,200,iter.max = 30)
print(all(sigInfo$sig_id==rownames(emb)))
sigInfo <- sigInfo %>% mutate(clusters=km$cluster)
df_info <- sigInfo %>% dplyr::select(sig_id,cmap_name,duplIdentifier,clusters) %>% unique()
saveRDS(df_info,'../results/clustering200_res.rds')

pointsInClusters <- df_info %>% group_by(clusters) %>% summarize(no_points=n_distinct(sig_id)) %>% ungroup()
hist(pointsInClusters$no_points,breaks=50)
df_info <- left_join(df_info,pointsInClusters)
stat3 <- df_info %>% filter(cmap_name=='STAT3')
stat3 <- stat3 %>% group_by(clusters) %>% mutate(proportion = n_distinct(sig_id)/nrow(stat3)) %>% ungroup()

uniqueClusterPopulations <- unique(stat3$no_points)

## Sample for each population 20k times and build a NULL distribution
NullProportions <- function(no_points,dfInfo,total_stat3=17,iters=20000){
  points <- rep(no_points,iters)
  prop <- NULL
  for (i in 1:iters){
    df_sample <-  sample_n(dfInfo, no_points)
    prop[i] <- nrow(df_sample %>% filter(cmap_name=='STAT3'))/total_stat3

  }
  return(data.frame(points,prop))
} 

library(doRNG)
df_nulls <- NULL

### calculate distances
df_nulls <- foreach(no = uniqueClusterPopulations) %dorng% {
  NullProportions(no_points = no ,dfInfo = df_info)
}
df<- do.call(rbind,df_nulls)

p_vals <- NULL
for (i in 1:length(uniqueClusterPopulations)){
  tmp_stat <- stat3 %>% filter(no_points==uniqueClusterPopulations[i])
  tmp_stat <- unique(tmp_stat$proportion)
  tmp <- df_nulls[[i]]
  p_vals[i] <- sum(tmp$prop>=tmp_stat)/nrow(tmp)
}
hist(p_vals,20)
p.adj <- p.adjust(p_vals,"bonferroni")
hist(p.adj,20)
print(uniqueClusterPopulations[which(p_vals<0.01)])
noPoints <- uniqueClusterPopulations[which(p_vals<0.01)]
clusterSTAT3 <- unique(stat3$clusters[which(stat3$no_points %in% noPoints)])
#clusterSTAT3 <- c(9)
stat3 <- left_join(stat3,data.frame(no_points=uniqueClusterPopulations,p_vals,p.adj))
saveRDS(stat3,'../results/stat3_latent_clusters_withpvals.rds')

# # Keep those that are similar with STAT3 in at least half of the cell-lines
# dist_filtered <- dist_filtered %>% filter(proportion>=0.5)
# latentgenes <- unique(c(dist_filtered$cmap_name.x,dist_filtered$cmap_name.y))
# latentgenes <- latentgenes[-which(latentgenes=='STAT3')]
# latentgenes <- data.frame(genes=latentgenes)
df_filtered <- df_info %>% filter(clusters %in% clusterSTAT3)
latentgenes <- unique(df_filtered$cmap_name)
latentgenes <- latentgenes[-which(latentgenes=='STAT3')]
latentgenes <- data.frame(genes=latentgenes)
geneInfo <- read.delim('../data/geneinfo_beta.txt')
latentgenes <- left_join(latentgenes,geneInfo,by=c("genes"="gene_symbol"))
df_gene_canditates <- readRDS('../results/gene_candidates.rds')
latentgenes <- latentgenes %>% filter(genes %in% df_gene_canditates$gene_symbol)
saveRDS(latentgenes,'../results/latentgenes.rds')

### Get consensus gene----
latentgenes <- readRDS('../results/latentgenes.rds') %>% dplyr::select(genes)
GSEAgenes <- readRDS('../results/GSEAgenes.rds')

consensus <- left_join(latentgenes,GSEAgenes)
consensus <- consensus %>% filter(!is.na(gene_id))

### For each GSEA candidate find out how many of the genes belonging ----
# in the path from source to STAT3 are in latent space important genes.

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
write_delim(interactions,'../results/releventInteractions.txt',delim = '\t')
write_delim(nodeAttr,'../results/releventNodeAttr.txt',delim = '\t')

# Find now latents in the gene paths to prioritize some candidates from GSEA
paths <- paths %>% filter(pathGenes!='STAT3')
latentGenesInPaths <- left_join(latentgenes,paths,by=c('genes'='pathGenes'))
latentGenesInPaths <- latentGenesInPaths %>% filter(!is.na(paths))
xlsx::write.xlsx(latentGenesInPaths,'../results/latentGenesInPaths.xlsx')
