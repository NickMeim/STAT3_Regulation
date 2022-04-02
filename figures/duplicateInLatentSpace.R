### Create all distance distribution plots from the embedded samples
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

# Load embeddings
emb <- read.csv('../../shrna_embs1024.csv')
emb <- emb %>% column_to_rownames('X')
emb <- as.matrix(emb)
gc()

# Keep some duplicates to be computationally feasible
duplicateSigs <- sigInfo %>% filter(dupl_counts>=5)
duplicatesIndentity <- unique(duplicateSigs$duplIdentifier)

emb <- emb[as.character(unique(duplicateSigs$sig_id)),]
gc()

# Euclidean distances or cosine similarities in latent space
#distance <- as.matrix(dist(emb, method = "euclidean",diag = F,upper = F))
library(lsa)
X <- t(emb)
distance <- cosine(X)
colnames(distance) <- rownames(emb)
rownames(distance) <- rownames(emb)

### Convert matrix into data frame
# Keep only unique (non-self) pairs
distance[lower.tri(distance,diag = T)] <- -100
dist <- reshape2::melt(distance)
dist <- dist %>% filter(value != -100)

# Keep only specific columns for the original data frame
duplicateSigs <- duplicateSigs %>% 
  dplyr::select(sig_id,cmap_name,tas,nsample,duplIdentifier)

# Merge meta-data info and distances values
dist <- left_join(dist,duplicateSigs,by = c("Var1"="sig_id"))
dist <- left_join(dist,duplicateSigs,by = c("Var2"="sig_id"))
dist <- dist %>% mutate(is_duplicate = (duplIdentifier.x==duplIdentifier.y))
dist <- dist %>% filter(!is.na(value))
#saveRDS(dist,'../results/sim_cosine_latent_5dupls.rds')
gc()

# Visualize distributions
df_dist <- dist %>%  
  mutate(is_duplicate=ifelse(is_duplicate==T,
                             'Duplicate Signatures','Random Signatures')) %>%
  mutate(is_duplicate = factor(is_duplicate,
                               levels = c('Random Signatures',
                                          'Duplicate Signatures')))

ggplot(df_dist,aes(x=value,color=is_duplicate,fill=is_duplicate)) +
  geom_density(alpha=0.2) +
  labs(col = 'Type',fill='Type',title="Duplicates seperation in the latent space",x="Euclidean distance", y = "Density")+
  theme_classic() + theme(text = element_text(size=10))

# Thresholds of tas numbers to split dataset.
# TAS number is a metric given by the L1000 platform,
# which signifys the strength of the signal and
# the quality of the data. The higher it is
# the higher the quality of the data.

tas_thresholds_lower <- c(0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5)
tas_thresholds_upper <- c(0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55)

# Initialize empty lists to store plots
# and number of data after filtering 
# based on TAS number
plotList <- NULL
p_vals <- NULL
for (i in 1:length(tas_thresholds_upper)){
  
  # Get the distances for these TAS numbers
  df_dist <- dist %>%  
    filter((tas.x>=tas_thresholds_lower[i] & tas.y>=tas_thresholds_lower[i])) %>%
    mutate(is_duplicate=ifelse(is_duplicate==T,
                               'Duplicate Signatures','Random Signatures')) %>%
    mutate(is_duplicate = factor(is_duplicate,
                                 levels = c('Random Signatures',
                                            'Duplicate Signatures')))
  
  # Perform t-test
  random <- df_dist %>% filter(is_duplicate=='Random Signatures')
  dupl <- df_dist %>% filter(is_duplicate!='Random Signatures')
  #p_vals[i] <- t.test(random$value,dupl$value,'greater')$p.value
  p_vals[i] <- t.test(random$value,dupl$value,'less')$p.value
  
  # Plot the distributions and store them forn now in the list using ggplot
  plotList[[i]] <- ggplot(df_dist,aes(x=value,color=is_duplicate,fill=is_duplicate)) +
    geom_density(alpha=0.2) +
    labs(col = 'Type',fill='Type',title="",x="Cosine Similarity", y = "Density")+
    theme_classic() + theme(text = element_text(size=10))
}
png(file="duplicate_vs_random_latent_cosine_distribution_shrna.png",width=16,height=8,units = "in",res=300)
p <- ggarrange(plotlist=plotList,ncol=4,nrow=2,common.legend = TRUE,legend = 'right',
               labels = c('TAS>=0.15','TAΣ>=0.2','TAS>=0.25','TAS>=0.3','TAS>=0.35',
                                 'TAS>=0.4','TAS>=0.45','TAS>=0.5'),
               font.label = list(size = 10, color = "black", face = "plain", family = NULL),
               hjust=-0.15)

annotate_figure(p, top = text_grob("Distributions of latent cosine similarities for different TAS thesholds", 
                                   color = "black",face = 'plain', size = 14))
dev.off()
# Adjust p-values with Bonferroni and plot them
p.adj <- p.adjust(p_vals,"bonferroni")
hist(p.adj)


