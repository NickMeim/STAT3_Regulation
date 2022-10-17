### Duplicate pathways distance distribution
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
sigInfo <- sigInfo %>% filter(tas>=0.3)

# Create identifier to signify duplicate
# signatures: meaning same drug, same dose,
# same time duration, same cell-type
sigInfo <- sigInfo %>%
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))

sigInfo <- sigInfo %>% group_by(duplIdentifier) %>% 
  mutate(dupl_counts = n()) %>% ungroup()

# Read distance pre-processed data (see preprocessingDistances.R 
# in preprocessing folderfor details)
dist <- readRDS('../results/dupl_shRNA_mean_dist_q1replicates_5dupls_long.rds')
dist <- dist %>% filter(is_duplicate==T)
dist <- dist %>% filter((Var1 %in% sigInfo$sig_id) & (Var2 %in% sigInfo$sig_id))
dist$value <- dist$value/2

# Plot the distributions and store them forn now in the list using ggplot
p <- ggplot(dist,aes(x=value)) +
    geom_density(alpha=0.2,color='blue',fill='blue') + 
    geom_vline(xintercept=0.3, linetype="dashed", color = "red")+
    labs(title="KEGG pathways GSEA distance of duplicate shRNAs",x="GSEA Distance", y = "Density")+
    xlim(c(0,1))+#xlim(c(min(df_dist$value),max(df_dist$value)))+
    theme_classic() + theme(text = element_text(size=20),
                            legend.position="none",plot.title=element_text(hjust = 0.5))
print(p)

# Save all subplots/distributions into one common figure
png(file="duplicate_distribution_kegg_shrna_tasfiltered.png",width=12,height=9,units = "in",res=300)
print(p)
dev.off()

# Cumulative distribution plot
png(file="duplicate_Cumulative_distribution_kegg_shrna_tasfiltered.png",width=12,height=9,units = "in",res=300)
plot(ecdf(dist$value),main='Cumulative probability distribution',
     xlab='GSEA distance')
dev.off()
