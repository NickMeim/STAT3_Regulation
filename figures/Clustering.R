library(tidyverse)

library(doFuture)
#options(future.globals.maxSize= 891289600)
#options(future.globals.maxSize= 4831838208)
# parallel set number of workers
cores <- 16
registerDoFuture()
plan(multiprocess,workers = cores)

# Load and run many clusters----
emb <- read.csv('../../shrna_embs1024.csv')
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

emb <- emb[which(rownames(emb) %in% sigInfo$sig_id),]

#### Clustering analysis ###
ks <- c(1,5,10,15,20,25,50,75,100,150,200)
ks <- c(ks,seq(250,2000,by=250),seq(2500,4500,by=500))
ParallelKmeansCluster <- function(kcenters,embeddings){
  dunn <- NULL
  for (j in 1:5){
    km <- kmeans(embeddings,kcenters,iter.max = 30)
    dunn[j] <- (km$tot.withinss)/(km$tot.withinss+km$betweenss)
  }
  mu_dunn <- mean(dunn)
  std_dunn <- sd(dunn)
  return(data.frame(kcenters,mu_dunn,std_dunn))
}

library(doRNG)
df_km <- NULL
print('Begin running kmeans....')
### calculate distances
df_km <- foreach(k = ks) %dorng% {
  ParallelKmeansCluster(kcenters = k ,embeddings = emb)
}
df_km <- do.call(rbind,df_km)
#saveRDS(df_km,'df_kelbow.rds')

## Elbow plot----
library(ggplot2)
library(KneeArrower)
ks <- seq(0,10000,1000)
ks[1] <- 1
dunn <- df_km$mu_dunn
k <- df_km$kcenters
#ind <- which(df_km$k==100)
xy <- findCutoff(k,dunn,method='first')
x <- 250
y <- 0.7507536
png('kmeans_elbow_plot.png',width=12,height=8,units = "in",res=300)
ggplot(df_km,aes(kcenters,mu_dunn)) +geom_line(size=1) +ylim(c(0,1))+
  geom_point(size=1) +scale_x_continuous(breaks = ks) + 
  geom_errorbar(aes(ymin=mu_dunn-std_dunn, ymax=mu_dunn+std_dunn), width=.5,
                position=position_dodge(.9))+
  geom_ribbon(aes(ymin = mu_dunn-std_dunn, ymax = mu_dunn+std_dunn), fill = "grey70")+
  geom_segment(aes(x=x,xend =x,y=0, yend=y),size=0.75,color="red",
               arrow = arrow(length = unit(0.3, "cm"))) + 
  annotate("text",x=900,y=0.33,label="Elbow point",size=5)+
  xlab('k-number of clusters') + ylab('Mean Dunn Index') +ggtitle('Elbow plot for optimal number of clusters')+
  theme(text = element_text(size=20),legend.position = "none",plot.title = element_text(hjust = 0.5))
dev.off()

hist(df_km$std_dunn,breaks=20,main = 'Distribution of standard deviation of dunn index',xlab='Dunn standard deviation')
hist(df_km$mu_dunn,breaks=20,main = 'Distribution of mean of dunn index',xlab='Dunn mean value')

# Percentage of total STAT3 profiles in clusters
library(ggplot2)
library(ggsignif)
stat3_cluster <- readRDS('../results/stat3_latent_clusters_withpvals.rds')
stat3_cluster <- stat3_cluster %>% select(clusters,p_vals,proportion,no_points) %>% unique()
stat3_cluster$clusters <- factor(stat3_cluster$clusters,
                                 levels =stat3_cluster$clusters[order(stat3_cluster$clusters)])

p <- ggplot(stat3_cluster,aes(x=clusters,y=proportion,fill=p_vals))  + geom_bar(stat="identity") +
  labs(fill='p.value',title = 'Proportion of total STAT3 profiles in the cluster')+
  xlab('Cluster id') + ylab('Proportion') +
  geom_signif(
    y_position = c(0.18, 0.18,0.18), xmin = c(1, 4,5), xmax = c(1, 4,5),
    annotation = c("**", "**","**"), tip_length = 0)+
  theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
print(p)

png('cluster_proportion_barplot.png',width=12,height=8,units = "in",res=300)
print(p)
dev.off()
