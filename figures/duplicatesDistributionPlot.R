### Create all distance distribution plots 
### betwwen tas=0.15 and 0.5.
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
sigInfo <- sigInfo %>% filter(pert_type=='trt_cp')
sigInfo <- sigInfo %>% filter(quality_replicates==1)

# Create identifier to signify duplicate
# signatures: meaning same drug, same dose,
# same time duration, same cell-type
sigInfo <- sigInfo %>%
  mutate(duplIdentifier = paste0(cmap_name,"_",pert_idose,"_",pert_itime,"_",cell_iname))

sigInfo <- sigInfo %>% group_by(duplIdentifier) %>% 
  mutate(dupl_counts = n()) %>% ungroup()

# Read distance pre-processed data (see preprocessingDistances.R 
# in preprocessing folderfor details)
dist <- readRDS('../preprocessing/preprocessed_data/dupl_mean_dist.rds')


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
num_of_sigs <- NULL

for (i in 1:length(tas_thresholds_upper)){
  
  # Calculate number of data remaining
  num_of_sigs[i] <- nrow(sigInfo %>% filter(tas>=tas_thresholds_lower[i]))
  
  # Get the distances for these TAS numbers
  df_dist <- dist %>%  
    filter((tas.x>=tas_thresholds_lower[i] & tas.y>=tas_thresholds_lower[i])) %>%
    mutate(is_duplicate=ifelse(is_duplicate==T,
                               'Duplicate Signatures','Random Signatures')) %>%
    mutate(is_duplicate = factor(is_duplicate,
                                 levels = c('Random Signatures',
                                            'Duplicate Signatures')))
  # The distance metric ranges from 0-2
  # Normalize it between 0-1
  df_dist$value <- df_dist$value/2
  
  
  # Plot the distributions and store them forn now in the list using ggplot
  plotList[[i]] <- ggplot(df_dist,aes(x=value,color=is_duplicate,fill=is_duplicate)) +
    geom_density(alpha=0.2) +
    labs(col = 'Type',fill='Type',title="",x="GSEA Distance", y = "Density")+
    xlim(c(0,1))+#xlim(c(min(df_dist$value),max(df_dist$value)))+
    theme_classic() + theme(text = element_text(size=10))
}

# Save all subplots/distributions into one common figure
png(file="duplicate_vs_random_distribution.png",width=16,height=8,units = "in",res=300)

p <- ggarrange(plotlist=plotList,ncol=4,nrow=2,common.legend = TRUE,legend = 'right',
               labels = paste0(c('TAS>=0.15','TAΣ>=0.2','TAS>=0.25','TAS>=0.3','TAS>=0.35',
                                 'TAS>=0.4','TAS>=0.45','TAS>=0.5'),',Number of signatures:',num_of_sigs),
               font.label = list(size = 10, color = "black", face = "plain", family = NULL),
               hjust=-0.15)

annotate_figure(p, top = text_grob("Distributions of GSEA distances for different TAS thesholds", 
                                   color = "black",face = 'plain', size = 14))
dev.off()

# Save the line-plot of how the number of availabe data
# changes with varying TAS thresholds

#Dummy dataframe to plot
df = data.frame(tas_thresholds_lower,num_of_sigs)

# Plot the decrease of aavailable data
# as we become more strict in their quality
png(file="numberOfSignatures_vs_TASthresholf.png",width=9,height=6,units = "in",res=300)
ggplot(df,aes(x=tas_thresholds_lower,y=num_of_sigs)) +
  labs(x="TAS number threshold", y = "Number of data points")+
  geom_point() +geom_line(linetype='dashed')
dev.off()
