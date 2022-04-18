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
pca.samples <- prcomp(as.matrix(emb[,2:(ncol(emb)-1)]),scale=F)

#Visualize explained variance from top 20 PCs
png('explained_variance_pcs_latent.png',width=9,height=6,units = "in",res=300)
fviz_eig(pca.samples,ncp=50)
dev.off()

#Built data frame with PCs, labels etc
df_pca <- pca.samples$x[,1:2]
df_pca <- data.frame(df_pca)
df_pca$shRNA <- factor(emb$shRNA,levels = c('other','STAT3'))

#PCA plot
pca_plot <- ggplot(df_pca,aes(x=PC1,y=PC2))+ 
  ggtitle('Principal Component Analysis of latent space vectors')
pca_plot <- pca_plot+geom_point(aes(col =factor(shRNA),size=factor(shRNA)),alpha=0.5) +
  scale_color_manual(values = c("other"="#BFBFBF","STAT3"="#E46D25")) + 
  theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
print(pca_plot)

### tsne ####
library(Rtsne)
perpl = DescTools::RoundTo(sqrt(nrow(emb)), multiple = 5, FUN = round)
#Use the above formula to calculate perplexity (perpl). But if perplexity is too large for the number of data you have define manually
#perpl=2
init_dim = 20
iter = 1000
emb_size = ncol(emb)-2
tsne_all <- Rtsne(as.matrix(emb[,2:(ncol(emb)-1)]), 
                  dims = 2, perplexity=perpl, 
                  verbose=TRUE, max_iter = iter,initial_dims = init_dim,check_duplicates = F)
df_all <- data.frame(V1 = tsne_all$Y[,1], V2 =tsne_all$Y[,2])
df_all$shRNA <- factor(emb$shRNA,levels = c('other','STAT3'))

gtsne <- ggplot(df_all, aes(V1, V2))+
  geom_point(aes(col =factor(shRNA),size=factor(shRNA)))  + labs(title="t-SNE plot of latent space") + xlab('Dim 1') + ylab('Dim 2')+
  scale_color_manual(values = c("other"="#BFBFBF","STAT3"="#E46D25")) + theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
print(gtsne)



# Keep some duplicates to be computationally feasible
duplicateSigs <- sigInfo %>% filter(dupl_counts>=5)
duplicatesIndentity <- unique(duplicateSigs$duplIdentifier)

emb_dupls <- emb %>% filter(X %in% duplicateSigs$sig_id)
emb_dupl_stat3 <- emb %>% filter(X %in% stat3$sig_id)
emb_dupls <- emb_dupls %>% filter(!(X %in% stat3$sig_id))

ind <- which((emb$X %in% duplicateSigs$sig_id) | (emb$X %in% stat3$sig_id))
df <- df_all[ind,]

gtsne <- ggplot(df, aes(V1, V2))+
  geom_point(aes(col =factor(shRNA),size=factor(shRNA)))  + labs(title="t-SNE plot of latent space") + xlab('Dim 1') + ylab('Dim 2')+
  scale_color_manual(values = c("other"="#BFBFBF","STAT3"="#E46D25")) + theme(text = element_text(size=15),plot.title = element_text(hjust = 0.5))
print(gtsne)