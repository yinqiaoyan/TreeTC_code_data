##############################################################
###################### Figure S8 #############################
##############################################################

rm(list=ls())
load("Sim_data.RData")
load("Sim_result_iteration.RData")

# Please first install the following packages
library(ggplot2)
library(Rcpp)
library(tidyverse)
library(mcclust)
library(parallel)
library(mvtnorm)
library(aricode)
library(clevr)



##############################################################
############### (a) Scatter plot of data #####################
##############################################################

## True cluster memberships
trueClIds <- rep(c(rep(c(1, 4, 6, 7, 2, 3), each=10), 
                   rep(c(1, 5, 8, 9, 2, 3), each=10)), 5)

plot_data <- as.data.frame(testData)  # data in subtree
plot_data$clIds_plot <- trueClIds
p2 <- ggplot(plot_data, aes(x = V1, y=V2)) +
  geom_point(aes(color = factor(clIds_plot)))

## True centers
trueCenters = matrix(c(0, 0,-5, -5,5, 5,-5, 5,5, -5,-5, 7.5,-7.5, 5,5, -7.5,7.5, -5),
                     ncol = 2, byrow = T)

plot_df = as.data.frame(trueCenters)
colnames(plot_df) = c("x", "y")
finalCenters_df = plot_df

coord_list = list(
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[2],
             y2 = finalCenters_df$y[2]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[3],
             y2 = finalCenters_df$y[3]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[4],
             y2 = finalCenters_df$y[4]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[5],
             y2 = finalCenters_df$y[5]),
  
  data.frame(x1 = finalCenters_df$x[4],
             y1 = finalCenters_df$y[4],
             x2 = finalCenters_df$x[6],
             y2 = finalCenters_df$y[6]),
  
  data.frame(x1 = finalCenters_df$x[4],
             y1 = finalCenters_df$y[4],
             x2 = finalCenters_df$x[7],
             y2 = finalCenters_df$y[7]),
  
  data.frame(x1 = finalCenters_df$x[5],
             y1 = finalCenters_df$y[5],
             x2 = finalCenters_df$x[8],
             y2 = finalCenters_df$y[8]),
  
  data.frame(x1 = finalCenters_df$x[5],
             y1 = finalCenters_df$y[5],
             x2 = finalCenters_df$x[9],
             y2 = finalCenters_df$y[9])
)

p3 <- p2 + geom_point(aes(x=x, y=y), size=5, data = plot_df)

for (i in 1:length(coord_list)){
  p3 = p3 + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = coord_list[[i]],
                         size = 0.5, linetype=5)
}

p3 <- p3 + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-9, 9), ylim = c(-9, 9)) +
  scale_x_continuous(minor_breaks = c(-7.5, -5, 0, 5, 7.5), breaks = c(-7.5, -5, 0, 5, 7.5)) +
  scale_y_continuous(minor_breaks = c(-7.5, -5, 0, 5, 7.5), breaks = c(-7.5, -5, 0, 5, 7.5)) +
  theme(axis.title = element_blank(),
        axis.text = element_text(size = 25),
        legend.title = element_text(size = 35),
        legend.text = element_text(size = 30),
        legend.key.size = unit(3,"line")) +
  guides(color = guide_legend(override.aes = list(size = 6) ) )

ggsave("sim_truth.png", p3, width = 10.5, height = 8.1, dpi = 100)



##############################################################
########### (b) Clustering result of TreeTC ##################
##############################################################
colorSet = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F",
             "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3", "#A020F0")

######################
### Source group 1 ###
######################
treeInd <- groups_subj[1]
ggFinal <- ggFinal1
subFinalClIds <- subFinalClIds1
subFinalCenters <- subFinalCenters1
subData <- Reduce(rbind, DataOrigin[finalSubjAssignments == treeInd])

colorSubset = colorSet[c(1, 2, 3, 4, 6, 7)]

finalClIds_plot = rep(0, length(subFinalClIds))
for (i in 1:6){
  inddd = which(subFinalClIds == i)
  if (i == 1) finalClIds_plot[inddd] = 1
  if (i == 2) finalClIds_plot[inddd] = 4
  if (i == 3) finalClIds_plot[inddd] = 6
  if (i == 4) finalClIds_plot[inddd] = 7
  if (i == 5) finalClIds_plot[inddd] = 2
  if (i == 6) finalClIds_plot[inddd] = 3
}

plot_data <- as.data.frame(subData)  
plot_data$clIds_plot <- finalClIds_plot

finalCenters_df = as.data.frame(subFinalCenters)
colnames(finalCenters_df) = c("x", "y")

finalCenters_HasData = subFinalCenters[which(subFinalCenters[, 3] == 1), ]
plot_df = as.data.frame(finalCenters_HasData)
colnames(plot_df) = c("x", "y")

## Segments between center points
coord_list = list(
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[3],
             y2 = finalCenters_df$y[3]),
  
  data.frame(x1 = finalCenters_df$x[3],
             y1 = finalCenters_df$y[3],
             x2 = finalCenters_df$x[5],
             y2 = finalCenters_df$y[5]),
  
  data.frame(x1 = finalCenters_df$x[3],
             y1 = finalCenters_df$y[3],
             x2 = finalCenters_df$x[6],
             y2 = finalCenters_df$y[6]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[8],
             y2 = finalCenters_df$y[8]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[20],
             y2 = finalCenters_df$y[20])
)

p2 <- ggplot(plot_data, aes(x = V1, y=V2)) +
  geom_point(aes(color=factor(clIds_plot))) + scale_color_manual(values = colorSubset)

p3 <- p2 + geom_point(aes(x=x, y=y), size=5, data = plot_df)

for (i in 1:length(coord_list)){
  p3 <- p3 + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = coord_list[[i]], 
                          size = 0.5)
}

p3 <- p3 + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-9, 9), ylim = c(-9, 9)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),  
        legend.text = element_text(size = 35),
        legend.key.size = unit(4,"line")) + 
  guides(color = guide_legend(override.aes = list(size = 8) ) )

ggsave("sim_group1.png", p3, width = 10.5, height = 8.1, dpi = 100)


######################
### Source group 2 ###
######################
treeInd <- groups_subj[2]
ggFinal <- ggFinal2
subFinalClIds <- subFinalClIds2
subFinalCenters <- subFinalCenters2
subData <- Reduce(rbind, DataOrigin[finalSubjAssignments == treeInd])

colorSubset = colorSet[c(1, 2, 3, 5, 8, 9, 10)]

finalClIds_plot = rep(0, length(subFinalClIds))
for (i in 1:7){
  inddd = which(subFinalClIds == i)
  if (i == 1) finalClIds_plot[inddd] = 1
  if (i == 2) finalClIds_plot[inddd] = 2
  if (i == 3) finalClIds_plot[inddd] = 5
  if (i == 4) finalClIds_plot[inddd] = 9
  if (i == 5) finalClIds_plot[inddd] = 8
  if (i == 6) finalClIds_plot[inddd] = 10
  if (i == 7) finalClIds_plot[inddd] = 3
}

plot_data <- as.data.frame(subData)  
plot_data$clIds_plot <- finalClIds_plot

finalCenters_df = as.data.frame(subFinalCenters)
colnames(finalCenters_df) = c("x", "y")

finalCenters_HasData = subFinalCenters[which(subFinalCenters[, 3] == 1), ]
plot_df = as.data.frame(finalCenters_HasData)
colnames(plot_df) = c("x", "y")

## Segments between center points
coord_list = list(
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[8],
             y2 = finalCenters_df$y[8]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[11],
             y2 = finalCenters_df$y[11]),
  
  data.frame(x1 = finalCenters_df$x[11],
             y1 = finalCenters_df$y[11],
             x2 = finalCenters_df$x[14],
             y2 = finalCenters_df$y[14]),
  
  data.frame(x1 = finalCenters_df$x[11],
             y1 = finalCenters_df$y[11],
             x2 = finalCenters_df$x[16],
             y2 = finalCenters_df$y[16]),
  
  data.frame(x1 = finalCenters_df$x[16],
             y1 = finalCenters_df$y[16],
             x2 = finalCenters_df$x[18],
             y2 = finalCenters_df$y[18]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[20],
             y2 = finalCenters_df$y[20])
)

p2 <- ggplot(plot_data, aes(x = V1, y=V2)) +
  geom_point(aes(color=factor(clIds_plot))) + scale_color_manual(values= colorSubset)

p3 <- p2 + geom_point(aes(x=x, y=y), size=5, data = plot_df)

for (i in 1:length(coord_list)){
  p3 <- p3 + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = coord_list[[i]], 
                          size = 0.5)
}

p3 <- p3 + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-9, 9), ylim = c(-9, 9)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),  
        legend.text = element_text(size = 35),
        legend.key.size = unit(4,"line")) + 
  guides(color = guide_legend(override.aes = list(size = 8) ) )    

ggsave("sim_group2.png", p3, width = 10.5, height = 8.1, dpi = 100)



##############################################################
############### (c) Measures comparison ######################
##############################################################
add_legend = function (...) {
  opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 
                                                                0, 0, 0), new = TRUE)
  on.exit(par(opar))
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend(...)
}

## Compute observation-level clustering measures
trueClIds <- rep(c(rep(c(1, 4, 6, 7, 2, 3), each=10), 
                   rep(c(1, 5, 8, 9, 2, 3), each=10)), 5)

# methods: TreeTC, CAM, GMM, HC, k-centroids
clIds_method <- data.frame(clIds_treetc = finalClIds,
                           clIds_cam    = read.csv(file = "sim_cam_obs.csv")$c,
                           clIds_gmm    = read.csv(file = "sim_gmm.csv")$c,
                           clIds_hc     = read.csv(file = "sim_hc.csv")$c,
                           clIds_kc     = read.csv(file = "sim_kc.csv")$c)

# compute ARI, NMI, and v-measure
ari_vec = numeric(5)
nmi_vec = numeric(5)
vm_vec  = numeric(5)

for (i in 1:5) {
  ari_vec[i] = round(ARI(clIds_method[, i], trueClIds), 3)
  nmi_vec[i] = round(NMI(clIds_method[, i], trueClIds), 3)
  vm_vec[i]  = round(v_measure(trueClIds, clIds_method[, i]), 3)
}

# MCMC traces
# clIds_trace <- Reduce(rbind, Map(function(x){x[[numOfTrees+1]]}, clIds), c())
# tmp_m = matrix(nrow = nrow(clIds), ncol = 3)
# for (i in 1:nrow(clIds)) {
#   tmp_m[i, 1] = ARI(clIds_trace[i, ], trueClIds)
#   tmp_m[i, 2] = NMI(clIds_trace[i, ], trueClIds)
#   tmp_m[i, 3] = v_measure(trueClIds, clIds_trace[i, ])
# }
# 
# SaveData = as.data.frame(tmp_m)
# write.csv(SaveData, file = "sim_trace_measures.csv", row.names = F)
# [Note] Since this step is time-consuming, we save the matrix tmp_m in "sim_trace_measures.csv".

measure_treetc <- c(ari_vec[1], nmi_vec[1], vm_vec[1])
measure_cam    <- c(ari_vec[2], nmi_vec[2], vm_vec[2])
measure_gmm    <- c(ari_vec[3], nmi_vec[3], vm_vec[3])
measure_hc     <- c(ari_vec[4], nmi_vec[4], vm_vec[4])
measure_kc     <- c(ari_vec[5], nmi_vec[5], vm_vec[5])
measure_trace  <- read.csv(file = "sim_trace_measures.csv")


## plot
png("sim_measures.png", width = 1000, height = 790)

par(cex.lab=1, cex.axis=3.4, mgp = c(1, 1.8, 0))
boxplot(measure_trace, outline=F,
        ylim=c(0.3,0.95),
        cex.axis=3.4,
        border=c('gray60'), col='gray90',
        xaxt = "n")
axis(1, at=1:3, labels=c("ARI", "NMI", "v_measure"))
points(measure_treetc, pch=21,cex = 4, bg= 'red', col = "red")
pointData = data.frame(x=c(1, 2.05, 3.05), y = measure_gmm)
points(pointData, pch=15,cex = 3.5, bg= 'black')
pointData = data.frame(x=c(1, 1.95, 2.95), y = measure_cam)
points(pointData, pch=23,cex = 3.5, bg= 'deepskyblue2', col = "deepskyblue2")
points(measure_hc, pch=24,cex = 3.5, bg= 'black')
points(measure_kc, pch=25,cex = 3.5, bg= 'black')

colors1 <- c("gray90",'red', 'deepskyblue2', 'black', 'black', "black")
colors2 <- c("gray",'red', 'deepskyblue2', 'black', 'black', "black")
add_legend(x=0.3, y=0.2, legend=c("traces", 'TreeTC', 'CAM', "GMM",
                                  'k-centroids','HC'),
           pch=c(22,21,23,15,25,24), inset = c(0.1,0.13), col=colors1,
           pt.bg=colors2, 
           horiz=F, bty='n', cex=4)

dev.off()



##############################################################
############# (d) Clustering result of CAM ###################
##############################################################
colorSet = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F",
             "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3")
subj_clu <- read.csv(file = "sim_cam_sub.csv")$subj_clu
clu <- read.csv(file = "sim_cam_obs.csv")$c

######################
### Source group 1 ###
######################
colorSubset = colorSet[1:7]

subj_group = 1
subj_inds = which(subj_clu == subj_group)
csum_data = c(0, cumsum(rep(60, 10)))
subData = Reduce(rbind, Map(function(i){testData[(csum_data[i]+1):csum_data[i+1], ]}, subj_inds), c())
subClu = c()
for (i in subj_inds){
  subClu = c(subClu, clu[(csum_data[i]+1):csum_data[i+1]])
}

finalClIds_plot = rep(0, length(subClu))
for (i in c(1:6, 9)){
  inddd = which(subClu == i)
  if (i == 1) finalClIds_plot[inddd] = 1
  if (i == 2) finalClIds_plot[inddd] = 7
  if (i == 3) finalClIds_plot[inddd] = 4
  if (i == 4) finalClIds_plot[inddd] = 6
  if (i == 5) finalClIds_plot[inddd] = 2
  if (i == 6) finalClIds_plot[inddd] = 3
  if (i == 9) finalClIds_plot[inddd] = 5
}

plot_dat = data.frame(x = subData[,1], y = subData[,2], group = finalClIds_plot)
p = ggplot(data = plot_dat, aes(x=x, y=y)) + 
  geom_point(aes(color=factor(group))) + scale_color_manual(values= colorSubset)

p <- p + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-9, 9), ylim = c(-9, 9)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),  
        legend.text = element_text(size = 35)) + 
  theme(legend.key.size = unit(4,"line")) +     
  guides(color = guide_legend(override.aes = list(size = 8) ) )  

ggsave("sim_cam_group1.png", p, width = 10.5, height = 8.1, dpi = 100)


######################
### Source group 2 ###
######################
colorSubset = colorSet[c(1:5, 8, 9)]

subj_group = 2
subj_inds = which(subj_clu == subj_group)
csum_data = c(0, cumsum(rep(60, 10)))
subData = Reduce(rbind, Map(function(i){testData[(csum_data[i]+1):csum_data[i+1], ]}, subj_inds), c())
subClu = c()
for (i in subj_inds){
  subClu = c(subClu, clu[(csum_data[i]+1):csum_data[i+1]])
}

finalClIds_plot = rep(0, length(subClu))
for (i in c(1, 3, 5:9)){
  inddd = which(subClu == i)
  if (i == 1) finalClIds_plot[inddd] = 1
  if (i == 3) finalClIds_plot[inddd] = 4
  if (i == 5) finalClIds_plot[inddd] = 2
  if (i == 6) finalClIds_plot[inddd] = 3
  if (i == 7) finalClIds_plot[inddd] = 9
  if (i == 8) finalClIds_plot[inddd] = 8
  if (i == 9) finalClIds_plot[inddd] = 5
}

plot_dat = data.frame(x = subData[,1], y = subData[,2], group = finalClIds_plot)
p = ggplot(data = plot_dat, aes(x=x, y=y)) + 
  geom_point(aes(color=factor(group))) + scale_color_manual(values= colorSubset)

p <- p + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-9, 9), ylim = c(-9, 9)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),  
        legend.text = element_text(size = 35)) + 
  theme(legend.key.size = unit(4,"line")) +     
  guides(color = guide_legend(override.aes = list(size = 8) ) )  

ggsave("sim_cam_group2.png", p, width = 10.5, height = 8.1, dpi = 100)



##############################################################
##################      (e) GMM      #########################
##############################################################
colorSet = c("#F8766D", "#D39200", "#93AA00", "#00BA38",
             "#00B9E3", "#619CFF", "#FF61C3", "#6A5ACD", "#B03060")

zz_ids <- read.csv(file = "sim_gmm.csv")$c
finalClIds_plot = rep(0, length(zz_ids))
for (i in 1:9){
  inddd = which(zz_ids == i)
  if (i == 1) finalClIds_plot[inddd] = "C1"
  if (i == 2) finalClIds_plot[inddd] = 4
  if (i == 3) finalClIds_plot[inddd] = 1
  if (i == 4) finalClIds_plot[inddd] = 7
  if (i == 5) finalClIds_plot[inddd] = 9
  if (i == 6) finalClIds_plot[inddd] = 2
  if (i == 7) finalClIds_plot[inddd] = "C2"
  if (i == 8) finalClIds_plot[inddd] = 6
  if (i == 9) finalClIds_plot[inddd] = 3
}

plot_dat = data.frame(x = testData[,1], y = testData[,2], group = finalClIds_plot)
p2 <- ggplot(plot_dat, aes(x = x, y = y)) +
  geom_point(aes(color=factor(group))) + scale_color_manual(values = colorSet)

p3 <- p2 + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-9, 9), ylim = c(-9, 9)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 35),
        legend.key.size = unit(3,"line")) +
  guides(color = guide_legend(override.aes = list(size = 7.5) ) )

ggsave("sim_gmm.png", p3, width = 10.5, height = 8.1, dpi = 100)



##############################################################
######### (f) Clustering result of baseline methods ##########
##############################################################

######################
#####     HC     #####
######################

colorSubset = c("#9C661F", "#ED9121", "#C77CFF", "#9590FF")

hc_ids <- read.csv(file = "sim_hc.csv")$c
finalClIds_plot = rep(0, length(hc_ids))
for (i in 1:5){
  inddd = which(hc_ids == i)
  if (i == 1) finalClIds_plot[inddd] = "C3"
  if (i == 2) finalClIds_plot[inddd] = "C4"
  if (i == 3) finalClIds_plot[inddd] = "C5"
  if (i == 4) finalClIds_plot[inddd] = "C6"
}

plot_data <- as.data.frame(testData) 
plot_data$clIds_plot <- finalClIds_plot 

p2 <- ggplot(plot_data, aes(x = V1, y=V2)) +
  geom_point(aes(color = factor(clIds_plot)))  + scale_color_manual(values= colorSubset)

p3 <- p2 + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-9, 9), ylim = c(-9, 9)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 35),  
        legend.text = element_text(size = 35),
        legend.key.size = unit(4,"line")) + 
  guides(color = guide_legend(override.aes = list(size = 8) ) )  

ggsave("sim_hc.png", p3, width = 10.5, height = 8.1, dpi = 100)


######################
#### k-centroids #####
######################

colorSubset = c("#F8766D", "#D39200", "#93AA00", "#C77CFF", "#9590FF")

kc_ids <- read.csv(file = "sim_kc.csv")$c
finalClIds_plot = rep(0, length(kc_ids))
for (i in 1:5){
  inddd = which(kc_ids == i)
  if (i == 1) finalClIds_plot[inddd] = "C6"
  if (i == 2) finalClIds_plot[inddd] = 1
  if (i == 3) finalClIds_plot[inddd] = 3
  if (i == 4) finalClIds_plot[inddd] = 2
  if (i == 5) finalClIds_plot[inddd] = "C5"
}

plot_data <- as.data.frame(testData) 
plot_data$clIds_plot <- finalClIds_plot 

p2 <- ggplot(plot_data, aes(x = V1, y=V2)) +
  geom_point(aes(color = factor(clIds_plot)))  + scale_color_manual(values= colorSubset)

p3 <- p2 + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-9, 9), ylim = c(-9, 9)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 35),  
        legend.text = element_text(size = 35),
        legend.key.size = unit(4,"line")) + 
  guides(color = guide_legend(override.aes = list(size = 8) ) )  

ggsave("sim_kc.png", p3, width = 10.5, height = 8.1, dpi = 100)


