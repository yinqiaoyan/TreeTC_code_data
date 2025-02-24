##############################################################
###################### Figure S9 #############################
##############################################################

rm(list=ls())
load("Real_result_iteration.RData")
rawData <- read.csv(file = "rawData.csv", header = F) 
rawData <- as.matrix(rawData)
indices  <- rawData[, dim(rawData)[2]]

numLen = 10  # show 10 images with the highest likelihood in one node
numSel = 20  # only select nodes with more than 20 images
csum_data = c(0, cumsum(numOfSubjData))


##############################################################
################### Image group 1 ############################
##############################################################
flag = 0
number = 0
topiidMat = c()
flagRecord = c()

treeInd <- groups_subj[1]
subFinalClIds <- subFinalClIds1
final_center <- subFinalCenters1
indd <- which(finalSubjAssignments == treeInd)
subData <- Reduce(rbind, Map(function(i){testData[(csum_data[i]+1):csum_data[i+1], ]}, indd), c())
indices_sub <- unlist(Map(function(i){indices[(csum_data[i]+1):csum_data[i+1]]}, indd))

act_center = final_center[final_center[, ncol(final_center)] == 1, -ncol(final_center)]
act_variance = finalSigma
nodeNames = rownames(act_center)
depths = (nchar(nodeNames) - 1) / 2

while (flag < max(subFinalClIds)) {
  ids = which(subFinalClIds == flag+1)
  
  if (length(ids) >= numSel) {
    flagRecord = c(flagRecord, length(ids))
    number = number + 1
    llhRecord <- rep(0, length(ids))
    for (i in 1:length(ids)) {
      llhRecord[i] <- sum(dnorm(testData[ids[i], ], 
                                mean = act_center[flag+1, ], 
                                sd = sqrt(act_variance * etaNormal^(depths[flag+1])),
                                log = TRUE))
    }
    
    top10Llh <- sort(llhRecord, decreasing = TRUE)[1:numLen]
    top10Ids <- rep(0, numLen)
    
    for (i in 1:numLen) {
      top10Ids[i] <- which(llhRecord == top10Llh[i])
    }
    
    top10Ids <- ids[top10Ids]
    top10iid <- indices_sub[top10Ids]
    topiidMat <- rbind(topiidMat, top10iid)
  }
  flag = flag + 1
}

write.table(topiidMat, file = "./TopClIds/group1.csv", 
            row.names = F, col.names = F, sep = ",")


## Use the following codes to plot the tree structure of the images
dataCounts = numeric(max(subFinalClIds))
for (i in 1:max(subFinalClIds)) {
  ids = which(subFinalClIds == i)
  dataCounts[i] = length(ids)
}
dataFrame = data.frame(nodeNames = nodeNames, counts = dataCounts)
dataFrame_20 = dataFrame[dataCounts >= numSel, ]
# dataFrame_20 saves the node indices where more than numSel = 20 
# data points locate.



##############################################################
################### Image group 2 ############################
##############################################################
flag = 0
number = 0
topiidMat = c()
flagRecord = c()

treeInd <- groups_subj[2]
subFinalClIds <- subFinalClIds2
final_center <- subFinalCenters2
indd = which(finalSubjAssignments == treeInd)
subData <- Reduce(rbind, Map(function(i){testData[(csum_data[i]+1):csum_data[i+1], ]}, indd), c())
indices_sub <- unlist(Map(function(i){indices[(csum_data[i]+1):csum_data[i+1]]}, indd))

act_center = final_center[final_center[, ncol(final_center)] == 1, -ncol(final_center)]
act_variance = finalSigma
nodeNames = rownames(act_center)
depths = (nchar(nodeNames) - 1) / 2

while (flag < max(subFinalClIds)) {
  ids = which(subFinalClIds == flag+1)
  
  if (length(ids) >= numSel) {
    flagRecord = c(flagRecord, length(ids))
    number = number + 1
    llhRecord <- rep(0, length(ids))
    for (i in 1:length(ids)) {
      llhRecord[i] <- sum(dnorm(testData[ids[i], ], 
                                mean = act_center[flag+1, ], 
                                sd = sqrt(act_variance * etaNormal^(depths[flag+1])),
                                log = TRUE))
    }
    
    top10Llh <- sort(llhRecord, decreasing = TRUE)[1:numLen]
    top10Ids <- rep(0, numLen)
    
    for (i in 1:numLen) {
      top10Ids[i] <- which(llhRecord == top10Llh[i])
    }
    
    top10Ids <- ids[top10Ids]
    top10iid <- indices_sub[top10Ids]
    topiidMat <- rbind(topiidMat, top10iid)
  }
  flag = flag + 1
}

write.table(topiidMat, file = "./TopClIds/group2.csv", 
            row.names = F, col.names = F, sep = ",")


## Use the following codes to plot the tree structure of the images
dataCounts = numeric(max(subFinalClIds))
for (i in 1:max(subFinalClIds)) {
  ids = which(subFinalClIds == i)
  dataCounts[i] = length(ids)
}
dataFrame = data.frame(nodeNames = nodeNames, counts = dataCounts)
dataFrame_20 = dataFrame[dataCounts >= numSel, ]
# dataFrame_20 saves the node indices where more than numSel = 20 
# data points locate.



##############################################################
################### Image group 3 ############################
##############################################################
flag = 0
number = 0
topiidMat = c()
flagRecord = c()

treeInd <- groups_subj[3]
subFinalClIds <- subFinalClIds3
final_center <- subFinalCenters3
indd = which(finalSubjAssignments == treeInd)
subData <- Reduce(rbind, Map(function(i){testData[(csum_data[i]+1):csum_data[i+1], ]}, indd), c())
indices_sub <- unlist(Map(function(i){indices[(csum_data[i]+1):csum_data[i+1]]}, indd))

act_center = final_center[final_center[, ncol(final_center)] == 1, -ncol(final_center)]
act_variance = finalSigma
nodeNames = rownames(act_center)
depths = (nchar(nodeNames) - 1) / 2

while (flag < max(subFinalClIds)) {
  ids = which(subFinalClIds == flag+1)
  
  if (length(ids) >= numSel) {
    flagRecord = c(flagRecord, length(ids))
    number = number + 1
    llhRecord <- rep(0, length(ids))
    for (i in 1:length(ids)) {
      llhRecord[i] <- sum(dnorm(testData[ids[i], ], 
                                mean = act_center[flag+1, ], 
                                sd = sqrt(act_variance * etaNormal^(depths[flag+1])),
                                log = TRUE))
    }
    
    top10Llh <- sort(llhRecord, decreasing = TRUE)[1:numLen]
    top10Ids <- rep(0, numLen)
    
    for (i in 1:numLen) {
      top10Ids[i] <- which(llhRecord == top10Llh[i])
    }
    
    top10Ids <- ids[top10Ids]
    top10iid <- indices_sub[top10Ids]
    topiidMat <- rbind(topiidMat, top10iid)
  }
  flag = flag + 1
}

write.table(topiidMat, file = "./TopClIds/group3.csv", 
            row.names = F, col.names = F, sep = ",")


## Use the following codes to plot the tree structure of the images
dataCounts = numeric(max(subFinalClIds))
for (i in 1:max(subFinalClIds)) {
  ids = which(subFinalClIds == i)
  dataCounts[i] = length(ids)
}
dataFrame = data.frame(nodeNames = nodeNames, counts = dataCounts)
dataFrame_20 = dataFrame[dataCounts >= numSel, ]
# dataFrame_20 saves the node indices where more than numSel = 20 
# data points locate.





