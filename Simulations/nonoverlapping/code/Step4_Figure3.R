##############################################################
###################### Figure 3 ##############################
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
trueCenters = matrix(c(0, 0,-12, -12,12, 12,-12, 12,12, -12,-12, 16,-16, 12,12, -16,16, -12),
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
  coord_cartesian(xlim=c(-18, 18), ylim = c(-18, 18)) +
  scale_x_continuous(minor_breaks = c(-16, -12, 0, 12, 16), breaks = c(-16, -12, 0, 12, 16)) +
  scale_y_continuous(minor_breaks = c(-16, -12, 0, 12, 16), breaks = c(-16, -12, 0, 12, 16)) +
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
             "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3")

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
  if (i == 2) finalClIds_plot[inddd] = 3
  if (i == 3) finalClIds_plot[inddd] = 4
  if (i == 4) finalClIds_plot[inddd] = 6
  if (i == 5) finalClIds_plot[inddd] = 7
  if (i == 6) finalClIds_plot[inddd] = 2
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
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[6],
             y2 = finalCenters_df$y[6]),
  
  data.frame(x1 = finalCenters_df$x[6],
             y1 = finalCenters_df$y[6],
             x2 = finalCenters_df$x[9],
             y2 = finalCenters_df$y[9]),
  
  data.frame(x1 = finalCenters_df$x[6],
             y1 = finalCenters_df$y[6],
             x2 = finalCenters_df$x[12],
             y2 = finalCenters_df$y[12]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[14],
             y2 = finalCenters_df$y[14])
)

p2 <- ggplot(plot_data, aes(x = V1, y=V2)) +
  geom_point(aes(color=factor(clIds_plot))) + scale_color_manual(values = colorSubset)

p3 <- p2 + geom_point(aes(x=x, y=y), size=5, data = plot_df)

for (i in 1:length(coord_list)){
  p3 <- p3 + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = coord_list[[i]], 
                          size = 0.5)
}

p3 <- p3 + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-18, 18), ylim = c(-18, 18)) +
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

colorSubset = colorSet[c(1, 2, 3, 5, 8, 9)]

finalClIds_plot = rep(0, length(subFinalClIds))
for (i in 1:6){
  inddd = which(subFinalClIds == i)
  if (i == 1) finalClIds_plot[inddd] = 1
  if (i == 2) finalClIds_plot[inddd] = 3
  if (i == 3) finalClIds_plot[inddd] = 2
  if (i == 4) finalClIds_plot[inddd] = 5
  if (i == 5) finalClIds_plot[inddd] = 8
  if (i == 6) finalClIds_plot[inddd] = 9
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
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[14],
             y2 = finalCenters_df$y[14]),
  
  data.frame(x1 = finalCenters_df$x[1],
             y1 = finalCenters_df$y[1],
             x2 = finalCenters_df$x[17],
             y2 = finalCenters_df$y[17]),
  
  data.frame(x1 = finalCenters_df$x[17],
             y1 = finalCenters_df$y[17],
             x2 = finalCenters_df$x[20],
             y2 = finalCenters_df$y[20]),
  
  data.frame(x1 = finalCenters_df$x[17],
             y1 = finalCenters_df$y[17],
             x2 = finalCenters_df$x[22],
             y2 = finalCenters_df$y[22])
)

p2 <- ggplot(plot_data, aes(x = V1, y=V2)) +
  geom_point(aes(color=factor(clIds_plot))) + scale_color_manual(values= colorSubset)

p3 <- p2 + geom_point(aes(x=x, y=y), size=5, data = plot_df)

for (i in 1:length(coord_list)){
  p3 <- p3 + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = coord_list[[i]], 
                          size = 0.5)
}

p3 <- p3 + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-18, 18), ylim = c(-18, 18)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),  
        legend.text = element_text(size = 35),
        legend.key.size = unit(4,"line")) + 
  guides(color = guide_legend(override.aes = list(size = 8) ) )  

ggsave("sim_group2.png", p3, width = 10.5, height = 8.1, dpi = 100)



##############################################################
############# (c) Unnormalized posterior #####################
##############################################################
png("sim_llh.png", width = 1000, height = 790)

par(oma=c(1,1,1,1), mar=c(5,5,5,5))
plot(ll_UnPost[200:length(ll_UnPost)], type="l",
     xlab = "",
     ylab = "",
     cex.axis = 2.8,
     xaxt="n")  
mtext(expression(paste("MCMC samples")), side = 1, line = 5.2, cex=4)
axis(1, at = c(1, 4801, 9801, 14801, 19801, 24801, 29801),
     labels=c("200", "5000", "10000","15000", "20000","25000", "30000"),
     mgp=c(1, 1.5, 0), cex.axis = 2.8)

dev.off()



##############################################################
############# (d) Clustering result of CAM ###################
##############################################################
colorSet = c("#F8766D", "#D39200", "#93AA00", "#00BA38", "#00C19F",
             "#00B9E3", "#619CFF", "#DB72FB", "#FF61C3", "#6A5ACD")
subj_clu <- read.csv(file = "sim_cam_sub.csv")$subj_clu
clu <- read.csv(file = "sim_cam_obs.csv")$c

######################
### Source group 1 ###
######################
colorSubset = colorSet[c(1,2,3,4,6,7)]

subj_group = 1
subj_inds = which(subj_clu == subj_group)
csum_data = c(0, cumsum(rep(60, 10)))
subData = Reduce(rbind, Map(function(i){testData[(csum_data[i]+1):csum_data[i+1], ]}, subj_inds), c())
subClu = c()
for (i in subj_inds){
  subClu = c(subClu, clu[(csum_data[i]+1):csum_data[i+1]])
}

finalClIds_plot = rep(0, length(subClu))
for (i in 1:6){
  inddd = which(subClu == i)
  if (i == 1) finalClIds_plot[inddd] = 1
  if (i == 2) finalClIds_plot[inddd] = 4
  if (i == 3) finalClIds_plot[inddd] = 6
  if (i == 4) finalClIds_plot[inddd] = 7
  if (i == 5) finalClIds_plot[inddd] = 2
  if (i == 6) finalClIds_plot[inddd] = 3
}

plot_dat = data.frame(x = subData[,1], y = subData[,2], group = finalClIds_plot)
p = ggplot(data = plot_dat, aes(x=x, y=y)) + 
  geom_point(aes(color=factor(group))) + scale_color_manual(values= colorSubset)

p <- p + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-18, 18), ylim = c(-18, 18)) +
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
colorSubset = colorSet[c(1,2,3,8,10)]

subj_group = 2
subj_inds = which(subj_clu == subj_group)
csum_data = c(0, cumsum(rep(60, 10)))
subData = Reduce(rbind, Map(function(i){testData[(csum_data[i]+1):csum_data[i+1], ]}, subj_inds), c())
subClu = c()
for (i in subj_inds){
  subClu = c(subClu, clu[(csum_data[i]+1):csum_data[i+1]])
}

finalClIds_plot = rep(0, length(subClu))
for (i in c(1,5,6,7,8)){
  inddd = which(subClu == i)
  if (i == 1) finalClIds_plot[inddd] = 1
  if (i == 5) finalClIds_plot[inddd] = 2
  if (i == 6) finalClIds_plot[inddd] = 3
  if (i == 7) finalClIds_plot[inddd] = "C1"
  if (i == 8) finalClIds_plot[inddd] = 8
}

plot_dat = data.frame(x = subData[,1], y = subData[,2], group = finalClIds_plot)
p = ggplot(data = plot_dat, aes(x=x, y=y)) + 
  geom_point(aes(color=factor(group))) + scale_color_manual(values= colorSubset)

p <- p + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-18, 18), ylim = c(-18, 18)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),  
        legend.text = element_text(size = 35)) + 
  theme(legend.key.size = unit(4,"line")) +     
  guides(color = guide_legend(override.aes = list(size = 8) ) )  

ggsave("sim_cam_group2.png", p, width = 10.5, height = 8.1, dpi = 100)



##############################################################
############### (e) Measures comparison ######################
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
# NOTE: Since this step is time-consuming, we save the matrix tmp_m in "sim_trace_measures.csv".

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
        ylim=c(0.4,1),
        cex.axis=3.4,
        border=c('gray60'), col='gray90',
        xaxt = "n")
axis(1, at=1:3, labels=c("ARI", "NMI", "v_measure"))

points(measure_treetc, pch=21,cex = 4, bg= 'red', col = "red")
pointData = data.frame(x=c(1.05, 2.05, 3.05), y = measure_gmm)
points(pointData, pch=15,cex = 3.5, bg= 'black')
pointData = data.frame(x=c(0.95, 1.95, 2.95), y = measure_cam)
points(pointData, pch=23,cex = 3.5, bg= 'deepskyblue2', col = "deepskyblue2")
pointData = data.frame(x=c(1.05, 2.05, 3.05), y = measure_hc)
points(pointData, pch=24,cex = 3.5, bg= 'black')
pointData = data.frame(x=c(0.95, 1.95, 2.95), y = measure_kc)
points(pointData, pch=25,cex = 3.5, bg= 'black')

colors1 <- c("gray90",'red', 'deepskyblue2', 'black', 'black', "black")
colors2 <- c("gray",'red', 'deepskyblue2', 'black', 'black', "black")
add_legend(x=0.3, y=0.2, legend=c("traces", 'TreeTC', 'CAM', "GMM",
                                  'k-centroids','HC'),
           pch=c(22,21,23,15,25,24), inset = c(0.1,0.13), col=colors1,
           pt.bg=colors2, 
           horiz=F, bty='n', cex=4)

dev.off()



##############################################################
##### (f) Clustering result of GMM and baseline methods ######
##############################################################

######################
######## GMM #########
######################

colorSet <- c("#F8766D", "#D39200", "#93AA00", "#00BA38",
              "#00B9E3", "#619CFF", "#DB72FB", "#6A5ACD", "#B03060")

zz_ids <- read.csv(file = "sim_gmm.csv")$c
finalClIds_plot = rep(0, length(zz_ids))
for (i in 1:9){
  inddd = which(zz_ids == i)
  if (i == 1) finalClIds_plot[inddd] = 7
  if (i == 2) finalClIds_plot[inddd] = 1
  if (i == 3) finalClIds_plot[inddd] = 3
  if (i == 4) finalClIds_plot[inddd] = 2
  if (i == 5) finalClIds_plot[inddd] = 4
  if (i == 6) finalClIds_plot[inddd] = "C2"
  if (i == 7) finalClIds_plot[inddd] = 6
  if (i == 8) finalClIds_plot[inddd] = 8
  if (i == 9) finalClIds_plot[inddd] = "C1"
}

plot_dat = data.frame(x = testData[,1], y = testData[,2], group = finalClIds_plot)
p2 <- ggplot(plot_dat, aes(x = x, y = y)) +
  geom_point(aes(color=factor(group))) + scale_color_manual(values = colorSet)

p3 <- p2 + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-18, 18), ylim = c(-18, 18)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),
        legend.text = element_text(size = 35),
        legend.key.size = unit(3,"line")) +
  guides(color = guide_legend(override.aes = list(size = 7.5) ) )

ggsave("sim_gmm.png", p3, width = 10.5, height = 8.1, dpi = 100)


######################
## HC / k-centroids ##
######################

# The result of hierarchical clustering and k-centroids are the same.
colorSubset = c("#F8766D", "#D39200", "#93AA00", "#9590FF", "#C77CFF")

hc_ids <- read.csv(file = "sim_hc.csv")$c
finalClIds_plot = rep(0, length(hc_ids))
for (i in 1:5){
  inddd = which(hc_ids == i)
  if (i == 1) finalClIds_plot[inddd] = 1
  if (i == 2) finalClIds_plot[inddd] = "C3"
  if (i == 3) finalClIds_plot[inddd] = 2
  if (i == 4) finalClIds_plot[inddd] = 3
  if (i == 5) finalClIds_plot[inddd] = "C4"
}

plot_data <- as.data.frame(testData) 
plot_data$clIds_plot <- finalClIds_plot 

p2 <- ggplot(plot_data, aes(x = V1, y=V2)) +
  geom_point(aes(color = factor(clIds_plot)))  + scale_color_manual(values= colorSubset)

p3 <- p2 + labs(color = "Clusters\n") +
  coord_cartesian(xlim=c(-18, 18), ylim = c(-18, 18)) +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_text(size = 40),  
        legend.text = element_text(size = 35)) + 
  theme(legend.key.size = unit(4,"line")) +      
  guides(color = guide_legend(override.aes = list(size = 8) ) )  

ggsave("sim_hckc.png", p3, width = 10.5, height = 8.1, dpi = 100)




