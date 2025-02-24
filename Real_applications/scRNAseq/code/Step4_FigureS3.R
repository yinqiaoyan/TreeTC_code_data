##############################################################
###################### Figure S3 #############################
##############################################################

rm(list=ls())
load("Real_data.RData")
load("Real_result_iteration.RData")

##############################################################
###################### Heatmaps ##############################
##############################################################
# Please first install the following packages
library(ggplot2)
library(tidyr)
library(tibble)
library(dplyr)
library(RColorBrewer)

uniqueCellTypes = sort(unique(annotation$cells))
maxCl = max(finalClIds)
n_clusters = maxCl
n_celltypes = length(uniqueCellTypes)

heat_data <- matrix(0, nrow=n_clusters, ncol=n_celltypes)
for (i in 1:n_clusters) {
  tmpids = which(finalClIds == i)
  tmpCT = dataFile$CellType[tmpids]
  heat_data[i, ] = as.numeric(unlist(Map(function(x){sum(tmpCT == x)}, uniqueCellTypes)))
}

## proportion
RowSums = rowSums(heat_data)
heat_data <- apply(heat_data, 2, function(x) {x / RowSums})

## ggplot for heatmap 
p <- heat_data %>%
  
  # Data wrangling
  as_tibble() %>%
  rowid_to_column(var="X") %>%
  gather(key="Y", value="Proportion", -1) %>%
  
  # Change Y to cell types
  mutate(Y=rep(uniqueCellTypes, each = maxCl)) %>%

  # Convert state to factor and reverse order of levels
  mutate(Y=factor(Y,levels=rev(sort(unique(Y))))) %>%

  # ggplot
  ggplot(aes(X, Y, fill= Proportion)) +
  geom_tile() +
  scale_fill_gradient(limits = c(0,1), low="white", high="#d53e4f") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 25),
        axis.text = element_text(size = 17),
        legend.title = element_text(size = 25),
        legend.text = element_text(size = 17),
        legend.key.size = unit(1.5,"line")) +
  labs(x = "Clusters", y = "Cell types") +
  scale_x_continuous(breaks=1:maxCl)

ggsave("Fig7.png", p, width = 12, height = 9, dpi = 100)

