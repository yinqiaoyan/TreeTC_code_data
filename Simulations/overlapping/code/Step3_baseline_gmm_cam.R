##############################################################
################ Methods for comparison ######################
##############################################################

rm(list=ls())
load("Sim_data.RData")

# Please first install the following packages
library(igraph)
library(ggplot2)
library(mcclust)
library(mclust)
library(cluster)
library(gplots)
library(ggvis)
library(Rcpp)


##############################################################
################### Baseline methods #########################
##############################################################

########################
###    Functions     ###
########################

##########################################################################
## These functions are developed by 
##
## Ke Yuan, Thomas Sakoparnig, Florian Markowetz and Niko Beerenwinkel. (2015) 
## BitPhylogeny: a probabilistic framework for reconstructing intra-tumor 
##    phylogenies. Genome Biology.
##########################################################################

## hierarchical clustering
get_label_hc <- function(x, K){
  
  dis <- dist(x, 'euclidean')
  hc_cand <- lapply(K, function(ii) cutree(hclust(dis), ii))
  hc_silhouette_res <- sapply(1:length(K),
                              function(ii)
                                summary(silhouette(hc_cand[[ii]], dis))$avg.width)
  idx <- which.max(hc_silhouette_res)
  hc_label <- hc_cand[[idx]]
  
  # compute the cluster mean and variance
  clone <- sapply(unique(hc_label), function(i) which(hc_label==i) )
  n <- length(clone)
  
  hc_sampleMean <- matrix(0, n, dim(x)[2])
  hc_sampleVar  <- matrix(0, n, dim(x)[2])  
  
  for (i in 1:n){
    idx <- clone[[i]]
    if (length(idx)==1){
      hc_sampleMean[i,] = as.matrix(x[idx,])
    }else{
      hc_sampleMean[i,] = colMeans(as.matrix(x[idx,]))
      hc_sampleVar [i,] = diag(cov(as.matrix(x[idx,])))
    }
  }
  return(list(label = hc_label, 
              center = hc_sampleMean, 
              hc_sampleVar = hc_sampleVar))
}


## k-centroids
get_label_kc <- function(x, K){
  
  dis <- dist(x, 'euclidean')
  
  kc_cand <- lapply(K, function(ii) pam( dis, ii) )
  
  kc_silhouette_res <- sapply(1:length(K),
                              function(ii)
                                summary( silhouette(kc_cand[[ii]]$clustering,
                                                    dis) )$avg.width )
  idx <- which.max( kc_silhouette_res )
  
  kc_label <- kc_cand[[idx]]$clustering
  kc_centroid <- x[kc_cand[[idx]]$medoids,]
  
  return(list(label = kc_label, center = kc_centroid))
}


## add legend on figure 
add_legend = function (...) {
  opar <- par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 
                                                                0, 0, 0), new = TRUE)
  on.exit(par(opar))
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend(...)
}

########################
### End of functions ###
########################

x <- testData
K <- 2:15

hcres <- get_label_hc(x, K)
kcres <- get_label_kc(x, K)

## Save results
SaveData = data.frame(x = testData[, 1], y = testData[, 2], c = hcres$label)
write.csv(SaveData, file = "sim_hc.csv", row.names = F)
SaveData = data.frame(x = testData[, 1], y = testData[, 2], c = kcres$label)
write.csv(SaveData, file = "sim_kc.csv", row.names = F)



##############################################################
######################     GMM    ############################
##############################################################

set.seed(30)
# Implement GMM five times under this seed and choose the best result

g <- 9
p <- ncol(testData)
n <- nrow(testData)
k_fit <- kmeans(testData, centers=g)

par <- vector("list", g)
par$pro <- k_fit$size/n
par$mean <- t(k_fit$centers)

variance <- mclustVariance("VII", d = p, G = g)
par$variance <- variance
par$variance$sigmasq <- c(4, 1, 1, 1, 1, 0.25, 0.25, 0.25, 0.25)

em_res <- em(modelName = "VII", data = testData, parameters = par)
zz <- em_res$z
zz_ids <- apply(zz, 1, function(x){which.max(x)})

## Save results
SaveData = data.frame(x = testData[, 1], y = testData[, 2], c = zz_ids)
write.csv(SaveData, file = "sim_gmm.csv", row.names = F)



##############################################################
######################     CAM    ############################
##############################################################

sourceCpp("../../CAM/newCAM.cpp")
source("../../CAM/newCAM.R")

Y <- testData
G <- rep(1:10, each=60)

set.seed(10)
# Implement CAM five times under this seed and choose the best result

ss = Sys.time()
OBJ <- CAM(y_obser = Y,
           y_group = G,
           K0 = 10,
           L0 = 10,
           prior = list(
             # hyperparameters NIG
             m0=mean(Y), k0=1/mean(var(Y)), a0=3, b0=1,
             # hyperparameters alpha and beta
             a_alpha=3, b_alpha = 3,
             a_beta =3, b_beta = 3 ),
           nsim = 30000,
           burn_in = 5000,
           thinning = 1,verbose = 1,
           fixedAB = T,
           kappa=0.5,cheap=T
)
ee = Sys.time()
exeTime = ee - ss

## Distributional Clusters
gt_distr   <- rep(1:2, 5)  # true unit group indices
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
maxGTD <- StatPerMeCo::Frobenius(
  matrix(0,nrow(DC_GT_PSM),ncol(DC_GT_PSM)),
  matrix(1,nrow(DC_GT_PSM),ncol(DC_GT_PSM)))

R <-  OBJ$Z_j
psm  <- PSM(R)
subj_clu  <- mcclust.ext::minVI(psm,method = "greedy")$cl
# https://github.com/sarawade/mcclust.ext
# devtools::install_github("sarawade/mcclust.ext")

## Observational Clusters
gt_distr   <- rep(c(rep(c(1, 4, 6, 7, 2, 3), each=10), 
                    rep(c(1, 5, 8, 9, 2, 3), each=10)), 5)  # true obs group indices
DC_GT_PSM  <- mcclust::comp.psm(rbind(gt_distr,gt_distr))
maxGTD <- StatPerMeCo::Frobenius(
  matrix(0,nrow(DC_GT_PSM),ncol(DC_GT_PSM)),
  matrix(1,nrow(DC_GT_PSM),ncol(DC_GT_PSM)))

R <-  OBJ$Csi_ij
psm  <- PSM(R)
clu  <- mcclust.ext::minVI(psm)$cl

## Save results
SaveData = data.frame(subj_clu)
write.csv(SaveData, file = "sim_cam_sub.csv", row.names = F)
SaveData = data.frame(x = testData[, 1], y = testData[, 2], c = clu)
write.csv(SaveData, file = "sim_cam_obs.csv", row.names = F)

