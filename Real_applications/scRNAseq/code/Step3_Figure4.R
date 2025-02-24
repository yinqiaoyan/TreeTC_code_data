##############################################################
###################### Figure 4 ##############################
##############################################################

rm(list=ls())
load("Real_result_iteration.RData")

# Please first install the following packages
library(ggplot2)

colorSet = c("#F8766D", "#C49A00", "#53B400", "#00C094", 
             "#00B6EB", "#A58AFF", "#FB61D7") 
finalClIds_all <- finalClIds
csum_data = c(0, cumsum(numOfSubjData))

##############################################################
################## Patient group 1 ###########################
##############################################################
# eight patients: 
# BMET1, BMET2, BMET3, BMET5, BMET6, BMET7, BMET8, BMET11.

treeInd <- groups_subj[1]
ggFinal <- ggFinal1
subFinalClIds <- subFinalClIds1
subFinalCenters <- subFinalCenters1

finalCenters_df = as.data.frame(subFinalCenters[, -ncol(subFinalCenters)])  # finalCenters_df used to draw line segment.
colnames(finalCenters_df) = c("x", "y")

finalCenters_HasData = subFinalCenters[c(1, 3:8), -ncol(subFinalCenters)]
plot_df = as.data.frame(finalCenters_HasData)  # plot_df used to draw the center points.
colnames(plot_df) = c("x", "y")

## Segments between center points
coord_list = list(
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[2],
             y2 = finalCenters_df$y[2]),
  
  data.frame(x1 = finalCenters_df$x[2],
             y1 = finalCenters_df$y[2],
             x2 = finalCenters_df$x[3],
             y2 = finalCenters_df$y[3]),
  
  data.frame(x1 = finalCenters_df$x[2],
             y1 = finalCenters_df$y[2],
             x2 = finalCenters_df$x[4],
             y2 = finalCenters_df$y[4]),
  
  data.frame(x1 = finalCenters_df$x[2],
             y1 = finalCenters_df$y[2],
             x2 = finalCenters_df$x[5],
             y2 = finalCenters_df$y[5]),
  
  data.frame(x1 = finalCenters_df$x[2],
             y1 = finalCenters_df$y[2],
             x2 = finalCenters_df$x[6],
             y2 = finalCenters_df$y[6]),
  
  data.frame(x1 = finalCenters_df$x[2],
             y1 = finalCenters_df$y[2],
             x2 = finalCenters_df$x[7],
             y2 = finalCenters_df$y[7]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[8],
             y2 = finalCenters_df$y[8])
)

## Draw figure for each patient in group 1.
for (indd in c(1, 3:9)){
  subData <- testData[(csum_data[indd]+1):csum_data[indd+1], ]

  plot_data <- as.data.frame(subData)
  plot_data$clIds <- finalClIds_all[(csum_data[indd]+1):csum_data[indd+1]]
  colorSubset <- colorSet
  
  p2 <- ggplot(plot_data, aes(x = V1, y=V2)) +
    geom_point(aes(color=factor(clIds))) + scale_color_manual(values= colorSubset)
  
  p3 = p2 + geom_point(aes(x=x, y=y), size=3, data = plot_df)
  
  for (i in 1:length(coord_list)){
    p3 = p3 + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = coord_list[[i]], 
                           size = 0.5)
  }
  
  plot_empty = data.frame(x = finalCenters_df$x[2], y = finalCenters_df$y[2])
  p3 = p3 + geom_point(aes(x, y), size=10, shape=1, data = plot_empty)
  
  p3 = p3 + labs(color = "Clusters\n") +
    coord_cartesian(xlim=c(-5.5, 10.3), ylim = c(-10.3, 8.8)) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_text(size = 40),
          legend.text = element_text(size = 35),
          legend.key.size = unit(4,"line")) +
    guides(color = guide_legend(override.aes = list(size = 8) ) )
  
  ggsave(paste0("bioData_subject", indd, ".png"), p3, width = 10.5, height = 8.1, dpi = 100)
}





##############################################################
################## Patient group 2 ###########################
##############################################################
# one patient: BMET10.

treeInd <- groups_subj[2]
ggFinal <- ggFinal2
subFinalClIds <- subFinalClIds2
subFinalCenters <- subFinalCenters2

finalCenters_df = as.data.frame(subFinalCenters[, -ncol(subFinalCenters)])  # finalCenters_df used to draw line segment.
colnames(finalCenters_df) = c("x", "y")

finalCenters_HasData = subFinalCenters[c(1,7,8), -ncol(subFinalCenters)]
plot_df = as.data.frame(finalCenters_HasData)  # plot_df used to draw the center points.
colnames(plot_df) = c("x", "y")

## Segments between center points
coord_list = list(
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[2],
             y2 = finalCenters_df$y[2]),
  
  data.frame(x1 = finalCenters_df$x[2],
             y1 = finalCenters_df$y[2],
             x2 = finalCenters_df$x[7],
             y2 = finalCenters_df$y[7]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[8],
             y2 = finalCenters_df$y[8])
)

## Draw figure for each patient in group 2.
for (indd in 2){
  subData <- testData[(csum_data[indd]+1):csum_data[indd+1], ]
  
  plot_data <- as.data.frame(subData)
  plot_data$clIds <- finalClIds_all[(csum_data[indd]+1):csum_data[indd+1]]
  colorSubset <- colorSet
  
  p2 <- ggplot(plot_data, aes(x = V1, y=V2)) +
    geom_point(aes(color=factor(clIds))) + scale_color_manual(values= colorSubset)
  
  p3 = p2 + geom_point(aes(x=x, y=y), size=3, data = plot_df)
  
  for (i in 1:length(coord_list)){
    p3 = p3 + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = coord_list[[i]], 
                           size = 0.5)
  }
  
  plot_empty = data.frame(x = finalCenters_df$x[2], y = finalCenters_df$y[2])
  p3 = p3 + geom_point(aes(x, y), size=10, shape=1, data = plot_empty)
  
  p3 = p3 + labs(color = "Clusters\n") +
    coord_cartesian(xlim=c(-5.5, 10.3), ylim = c(-10.3, 8.8)) +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.title = element_text(size = 40),
          legend.text = element_text(size = 35),
          legend.key.size = unit(4,"line")) +
    guides(color = guide_legend(override.aes = list(size = 8) ) )
  
  ggsave(paste0("bioData_subject", indd, ".png"), p3, width = 10.5, height = 8.1, dpi = 100)
}
