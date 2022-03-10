library(tidyverse)
library(cmapR)
library(org.Hs.eg.db)
library(rhdf5)
library(MLmetrics)

library(doFuture)

# parallel set number of workers
cores <- 15
registerDoFuture()
plan(multiprocess,workers = cores)

### Load data and keep only well-inferred and landmark genes----------------------------------------------------
geneInfo <- read.delim('../data/geneinfo_beta.txt')
geneInfo <-  geneInfo %>% filter(feature_space != "inferred")
# Keep only protein-coding genes
geneInfo <- geneInfo %>% filter(gene_type=="protein-coding")

#Load signature info and split data to high quality replicates and low quality replicates
sigInfo <- read.delim('../data/siginfo_beta.txt')
sigInfo <- sigInfo %>% mutate(quality_replicates = ifelse(is_exemplar_sig==1 & qc_pass==1 & nsample>=3,1,0))
sigInfo <- sigInfo %>% filter(pert_type=='trt_cp')
sigInfo <- sigInfo %>% filter(quality_replicates==1)

#Check correlation of duplicates CV and duplicates correlation with tas and nsamples
sigInfo <- sigInfo %>% mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))
sigInfo <- sigInfo %>% group_by(duplIdentifier) %>% mutate(dupl_counts = n()) %>% ungroup()

### Load gene expression data -----------------------------------------------------------------------------------

# Split sigs to run in parallel
sigIds <- unique(sigInfo$sig_id)
sigList <-  split(sigIds, ceiling(seq_along(sigIds)/ceiling(length(sigIds)/cores)))

# Parallelize parse_gctx
ds_path <- '../data/level5_beta_trt_cp_n720216x12328.gctx'

parse_gctx_parallel <- function(path ,rid,cid){
  gctx_file <- parse_gctx(path ,rid = rid,cid = cid)
  return(gctx_file@mat)
}
cmap_gctx <- foreach(sigs = sigList) %dopar% {
  parse_gctx_parallel(ds_path ,rid = unique(as.character(geneInfo$gene_id)),cid = sigs)
}
cmap <-do.call(cbind,cmap_gctx)
#saveRDS(cmap,'cmap_gex_compounds_qrep1_exemplar_rep3.rds')


### Check correlation of TAS with duplicates correlation and noise-------------------------------------------------
duplicateSigs <- sigInfo %>% filter(dupl_counts>1)
duplicatesIndentity <- unique(duplicateSigs$duplIdentifier)

# Calculate gsea distances 
library(doRNG)
# run distances
gex_dupls <- cmap[,unique(duplicateSigs$sig_id)]
thresholds <- c(30,50,100,200,300,400,500,600,700,800,900,1000)
dist_all_dupls <- NULL
### calculate distances
dist_all_dupls <- foreach(thres = thresholds) %dorng% {
  distance_scores(num_table = gex_dupls ,threshold_count = thres,names = colnames(gex_dupls))
}
distance <- do.call(cbind,dist_all_dupls)
distance <- array(distance,c(dim=dim(dist_all_dupls[[1]]),length(dist_all_dupls)))
mean_dist <- apply(distance, c(1,2), mean, na.rm = TRUE)
colnames(mean_dist) <- colnames(gex_dupls)
rownames(mean_dist) <- colnames(gex_dupls)
#saveRDS(mean_dist,'dupl_mean_dist.rds')
#mean_dist <- readRDS('dupl_mean_dist.rds')

### Compare distributions of duplicates vs random
mean_dist[lower.tri(mean_dist,diag = T)] <- -100
dist <- reshape2::melt(mean_dist)
dist <- dist %>% filter(value != -100)

duplicateSigs <- duplicateSigs %>% dplyr::select(sig_id,cmap_name,tas,nsample,duplIdentifier)

dist <- left_join(dist,duplicateSigs,by = c("Var1"="sig_id"))
dist <- left_join(dist,duplicateSigs,by = c("Var2"="sig_id"))
dist <- dist %>% mutate(is_duplicate = (duplIdentifier.x==duplIdentifier.y))
dist <- dist %>% mutate(AverageTAS = 0.5*(tas.x+tas.y))
dist <- dist %>% mutate(TasCV=sqrt((tas.x-AverageTAS)^2+(tas.y-AverageTAS)^2)/AverageTAS)

### Check all plot betwwen tas=0.15 and 0.5 ###
library(ggplot2)
library(ggpubr)
tas_thresholds_lower <- c(0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)
tas_thresholds_upper <- c(0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55)
differenceMeans <- NULL
tStatistic <- NULL
pValue <- NULL
plotList <- NULL
num_of_sigs <- NULL
for (i in 1:length(tas_thresholds_upper)){
  num_of_sigs[i] <- nrow(sigInfo %>% filter(tas>=tas_thresholds_lower[i]))
  dupl <- distDuplicate %>%
    filter((tas.x>=tas_thresholds_lower[i] & tas.y>=tas_thresholds_lower[i]))
  random <- distRandom %>% 
    filter((tas.x>=tas_thresholds_lower[i] & tas.y>=tas_thresholds_lower[i]))
  df_dist <- dist %>%  
    filter((tas.x>=tas_thresholds_lower[i] & tas.y>=tas_thresholds_lower[i])) %>%
    mutate(is_duplicate=ifelse(is_duplicate==T,'Duplicate Signatures','Random Signatures')) %>%
    mutate(is_duplicate = factor(is_duplicate,levels = c('Random Signatures','Duplicate Signatures')))
  df_dist$value <- df_dist$value/2
  test <- t.test(dupl$value/2,dist$value/2,alternative = 'less')
  tStatistic[i] <- test$statistic
  pValue[i] <- test$p.value
  differenceMeans[i] <- test$estimate['mean of y'] - test$estimate['mean of x']
  plotList[[i]] <- ggplot(df_dist,aes(x=value,color=is_duplicate,fill=is_duplicate)) +geom_density(alpha=0.2) +
    labs(col = 'Type',fill='Type',title="",x="GSEA Distance", y = "Density")+
    xlim(c(0,1))+#xlim(c(min(df_dist$value),max(df_dist$value)))+
    theme_classic() + theme(text = element_text(size=10))
}
png(file="duplicate_vs_random_distribution.png",width=16,height=8,units = "in",res=300)
p <- ggarrange(plotlist=plotList,ncol=4,nrow=2,common.legend = TRUE,legend = 'right',
          labels = paste0(c('TAS>=0.15','TAΣ>=0.2','TAS>=0.25','TAS>=0.3','TAS>=0.35',
                     'TAS>=0.4','TAS>=0.45','TAS>=0.5'),',Number of signatures:',num_of_sigs),
          font.label = list(size = 10, color = "black", face = "plain", family = NULL),
          hjust=-0.15)
annotate_figure(p, top = text_grob("Distributions of GSEA distances for different TAS thesholds", 
                                      color = "black",face = 'plain', size = 14))
dev.off()
df = data.frame(tas_thresholds_lower,num_of_sigs)
png(file="numberOfSignatures_vs_TASthresholf.png",width=9,height=6,units = "in",res=300)
ggplot(df,aes(x=tas_thresholds_lower,y=num_of_sigs)) +geom_point() +geom_line(linetype='dashed')
dev.off()