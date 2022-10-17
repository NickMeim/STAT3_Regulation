# Plot STAT3 in t-SNE space and then 1+9 random duplicates 
# Visualize the latent space in 2D
library(tidyverse)
library(ggplot2)
library(ggpubr)

#### Load meta-data of the L1000 platform
### It is the same as in pre-procesing

# Check L1000 documentation for information.

#Load signature info and split data to 
#high quality replicates and low quality replicates
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

stat3 <- sigInfo %>% filter(cmap_name=='STAT3')

# Load embeddings
emb <- read.csv('../../shrna_embs1024.csv')
emb <- emb %>% mutate(shRNA='other')
emb$shRNA[which(emb$X %in% stat3$sig_id)] <- 'STAT3'
#emb <- emb %>% column_to_rownames('X')
#emb <- as.matrix(emb)
gc()

# Visualize PCA and t-SNE STAT3 vs others

#pca analysis
library(factoextra)

#Run pca
pca.samples <- prcomp(as.matrix(emb[,2:(ncol(emb)-1)]),scale=T)

#Visualize explained variance from top 20 PCs
png('explained_variance_pcs_latent.png',width=9,height=6,units = "in",res=300)
fviz_eig(pca.samples,ncp=50)
dev.off()

#Built data frame with PCs, labels etc
df_pca <- pca.samples$x[,1:3]
df_pca <- data.frame(df_pca)
df_pca$shRNA <- factor(emb$shRNA,levels = c('other','STAT3'))

# library("gg3D")
# 
# #PCA plot
# pca_plot <- ggplot(df_pca,aes(x=PC1,y=PC2,z=PC3,col =shRNA))+ 
#   ggtitle('Principal Component Analysis of latent space vectors') +
#   scale_color_manual(values = c("other"="#BFBFBF","STAT3"="#E46D25"))+
#   theme_void() +
#   labs_3D(labs=c("PC1", "PC2", "PC3"),
#           angle=c(0,0,0),
#           hjust=c(0,2,2), 
#           vjust=c(2,2,-1))+
#   axes_3D() +
#   stat_3D()+
#   theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
# print(pca_plot)
# png('pca_3d_stat3.png',width=12,height=8,units = "in",res=300)
# print(pca_plot)
# dev.off()
pca_plot <- ggplot(df_pca, aes(PC1, PC2),alpha=0.2)+
  geom_point(aes(col =shRNA),alpha = 0.3)  + labs(title="PCA plot of latent space with STAT3 data-points") +
  xlab('PC1') + ylab('PC2')+
  scale_color_manual(values =c("other"="#BFBFBF","STAT3"="#E46D25")) +
  geom_point(data = subset(df_pca, shRNA != 'other'),aes(x = PC1, y = PC2, color = shRNA))+
  theme_minimal()+
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5))
print(pca_plot)
png('pca_2d_stat3.png',width=12,height=8,units = "in",res=300)
print(pca_plot)
dev.off()

df_cluster <- readRDS('../results/clustering200_res.rds')
df_pca$cluster <- df_cluster$clusters
releventClusters <- c(25,73,91)
df_pca <- df_pca %>% mutate(cluster=ifelse(cluster %in% releventClusters,cluster,'other'))
df_pca$cluster <- factor(df_pca$cluster, levels = c('25','73','91','other'))


# pca_plot <- ggplot(df_pca,aes(x=PC1,y=PC2,z=PC3,col =cluster))+ 
#   ggtitle('Principal Component Analysis of latent space vectors') +
#   scale_color_manual(values = c("25"="#E46D25","73"="#008080","91"="#800080","other"="#BFBFBF"))+
#   theme_void() +
#   labs_3D(labs=c("PC1", "PC2", "PC3"),
#           angle=c(0,0,0),
#           hjust=c(0,2,2), 
#           vjust=c(2,2,-1))+
#   axes_3D() +
#   stat_3D()+
#   theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
# 
# library(plotly)
# 
# fig <- plot_ly(df_pca, x = ~PC1, y = ~PC2, z = ~PC3,size=1,color=~shRNA,
#                colors = c("other"="#BFBFBF","STAT3"="#E46D25"),
#                mode   = 'markers',
#                type='scatter3d')
# fig %>% layout(title = "PCA 3D scatter plot of latent embeddings", 
#                scene = list(camera = list(eye = list(x = -1.3, y = 1.3, z = 1))))
# 
# df_pca <- df_pca %>% mutate(opacity=ifelse(cluster=='other',0.01,1))
# plot_ly(df_pca, x = ~PC1, y = ~PC2, z = ~PC3,size=1,color=~cluster,opacity = ~opacity,
#         
#         colors = c("25"="#E46D25","73"="#008080","91"="#800080","other"="#BFBFBF"),
#         mode   = 'markers',
#         type='scatter3d')
# 
# 
# 
# print(pca_plot)
# png('pca_3d_clusters.png',width=12,height=8,units = "in",res=300)
# print(pca_plot)
# dev.off()

pca_plot <- ggplot(df_pca, aes(PC1, PC2),alpha=0.2)+
  geom_point(aes(col =cluster),alpha = 0.3)  + labs(title="PCA plot of latent space with STAT3-relevant clusters") +
  xlab('PC1') + ylab('PC2')+
  scale_color_manual(values = c("25"="#E46D25","73"="#008080","91"="#800080","other"="#BFBFBF")) +
  geom_point(data = subset(df_pca, cluster != 'other'),aes(x = PC1, y = PC2, color = cluster))+
  theme_minimal()+
  theme(text = element_text(size=20),plot.title = element_text(hjust = 0.5))
print(pca_plot)
png('pca_2d_clusters.png',width=12,height=8,units = "in",res=300)
print(pca_plot)
dev.off()
### tsne ####
library(Rtsne)
#perpl = DescTools::RoundTo(sqrt(nrow(emb)), multiple = 5, FUN = round)
perpl= 250
init_dim = 10
iter = 2000
emb_size = ncol(emb)-2
set.seed(42)
tsne_all <- Rtsne(as.matrix(emb[,2:(ncol(emb)-1)]), 
                  dims = 2, perplexity=perpl, 
                  verbose=TRUE, 
                  max_iter = iter,
                  initial_dims = init_dim,
                  check_duplicates = F,
                  normalize = F,pca_scale = T,
                  num_threads = 15)
df_all <- data.frame(V1 = tsne_all$Y[,1], V2 =tsne_all$Y[,2])
df_all$shRNA <- factor(emb$shRNA,levels = c('other','STAT3'))

gtsne <- ggplot(df_all, aes(V1, V2),alpha=0.2)+
  geom_point(aes(col =factor(shRNA)))  + labs(title="t-SNE plot of latent space") + 
  xlab('Dim 1') + ylab('Dim 2')+
  scale_color_manual(values = c("other"="#BFBFBF","STAT3"="#E46D25")) + 
  geom_point(data = subset(df_all, shRNA == 'STAT3'),aes(x = V1, y = V2, color = factor(shRNA)))+
  theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
print(gtsne)

df_cluster <- readRDS('../results/clustering200_res.rds')
df_all$cluster <- df_cluster$clusters
releventClusters <- c(25,73,91)
df_all <- df_all %>% mutate(cluster=ifelse(cluster %in% releventClusters,cluster,'other'))
df_all$cluster <- factor(df_all$cluster, levels = c('25','73','91','other'))

gtsne <- ggplot(df_all, aes(V1, V2),alpha=0.2)+
  geom_point(aes(col =factor(cluster)),alpha = 0.3)  + labs(title="t-SNE plot of latent space with STAT3-relevant clusters") +
  xlab('Dim 1') + ylab('Dim 2')+
  scale_color_manual(values = c("25"="#E46D25","73"="#008080","91"="#800080","other"="#BFBFBF")) +
  geom_point(data = subset(df_all, cluster != 'other'),aes(x = V1, y = V2, color = cluster))+
  theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
print(gtsne)

### umap ###
library(umap)
library(plyr)
map <- umap(as.matrix(emb[,2:(ncol(emb)-1)]), n_components = 3)
df_map <- data.frame(V1 = map$layout[,1], V2 = map$layout[,2],V3= map$layout[,3])
df_map$shRNA <- factor(emb$shRNA,levels = c('STAT3','other'))
df_map$cluster <- df_cluster$clusters
releventClusters <- c(25,73,91)
df_map <- df_map %>% mutate(cluster=ifelse(cluster %in% releventClusters,cluster,'other'))
df_map$cluster <- factor(df_map$cluster, levels = c('25','73','91','other'))

gg_map <- ggplot(df_map, aes(V1, V2))+
  geom_point(aes(col =factor(cluster)),alpha=0.2)  + labs(title="UMAP plot of latent space with STAT3-relevant clusters") +
  xlab('Dim 1') + ylab('Dim 2')+ xlim(c(-6,6)) + ylim(c(-6,6))+
  scale_color_manual(values = c("25"="#E46D25","73"="#008080","91"="#800080","other"="#BFBFBF")) +
  geom_point(data = subset(df_map, cluster != 'other'),aes(x = V1, y = V2, color = cluster))+
  theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
print(gg_map)

gg_map <- ggplot(df_map, aes(V1, V2))+
  geom_point(aes(col =shRNA),alpha=0.2)  + labs(title="UMAP plot of latent space with STAT3-relevant clusters") +
  xlab('Dim 1') + ylab('Dim 2')+ xlim(c(-6,6)) + ylim(c(-6,6))+
  scale_color_manual(values = c("other"="#BFBFBF","STAT3"="#E46D25")) + 
  geom_point(data = subset(df_map, shRNA == 'STAT3'),aes(x = V1, y = V2, color = shRNA))+
  theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
print(gg_map)

# UMAP 3D
gg_map <- ggplot(df_map,aes(x=V1,y=V2,z=V3,col =cluster))+ 
  ggtitle('UMAP of latent space vectors') +
  scale_color_manual(values = c("25"="#E46D25","73"="#008080","91"="#800080","other"="#BFBFBF"))+
  theme_void() +
  labs_3D(labs=c("Dim 1", "Dim 2", "Dim 3"),
          angle=c(0,0,0),
          hjust=c(0,2,2), 
          vjust=c(2,2,-1))+
  axes_3D() +
  stat_3D()+
  theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
print(gg_map)