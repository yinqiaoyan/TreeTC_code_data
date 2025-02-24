##############################################################
############# Implement TreeTC and conduct MCMC ##############
##############################################################

rm(list=ls())

# Please first install the following packages
library(igraph)
library(ggplot2)

# Import functions -------------------------------------------------------------
source('../../TreeTC_class/util.R', echo=F)  
source('../../TreeTC_class/distn.R', echo=F)
source('../../TreeTC_class/mvnorm.R', echo=F)
source('../../TreeTC_class/auxiliary.R', echo=F)

source('../../TreeTC_class/node_normal.R', echo=F)
source('../../TreeTC_class/node.R', echo=F)
source('../../TreeTC_class/tssb_inference_mcmc.R', echo=F)
source('../../TreeTC_class/tssb.R', echo=F)


# Import input data ------------------------------------------------------------
rawData <- read.csv(file = "rawData.csv", header = F) 
rawData <- as.matrix(rawData)
testData <- rawData[, -dim(rawData)[2]]


# Initialize model and training ------------------------------------------------
fixSeedForTrain <- 111
set.seed(fixSeedForTrain)

priorTypeSigma = "invGamma"
priorTypeDrift = "unif"

priorSigmaScale = 0.1

priorDriftMin = 0.1
priorDriftMax = 1

numOfTrees = 15
maxDepthYuan = 10

etaNormal = 0.5
etaTheta  = 0.5

dpKappa = 0.5
dpAlpha0 = 1
dpLambda0 = 0.5
dpGamma0 = 0.5
dpAlpha1 = 0.5
dpGamma1 = 0.5

numOfSubjects = 30
numOfSubjData = rep(nrow(testData) / numOfSubjects, numOfSubjects)


q0 <- Normal$new(priorSigmaScale = priorSigmaScale,
                 priorSigmaMin = NULL,
                 priorSigmaMax = NULL,
                 priorDriftScale = NULL,
                 priorDriftMin = priorDriftMin,
                 priorDriftMax = priorDriftMax,
                 etaNormal = etaNormal,
                 etaTheta  = etaTheta,
                 dataDims = ncol(testData),
                 priorTypeSigma = priorTypeSigma,
                 priorTypeDrift = priorTypeDrift,
                 numOfTrees = numOfTrees)


tssb <- TssbMCMC$new(q0, data = testData,
                     numOfSubjects = numOfSubjects,
                     numOfSubjData = numOfSubjData,
                     dpKappa = dpKappa,
                     dpAlpha0 = dpAlpha0, dpGamma0 = dpGamma0, dpLambda0 = dpLambda0,
                     dpAlpha1 = dpAlpha1, dpGamma1 = dpGamma1,
                     maxDepth = maxDepthYuan,
                     initType = "unif")



# Post Update ------------------------------------------------------------------
numOfMCMC <- 10000
burnIn <- 5000


# Start Training ---------------------------------------------------------------
ll_UnPost = array(0, dim = numOfMCMC)
gg = vector(mode = "list", length = numOfTrees + 1)
gg <- Reduce(rbind, Map(function(i) {list(gg)}, seq_len(numOfMCMC - burnIn)))

clIds <- vector(mode = "list", length = numOfTrees+1)
clIds <- Reduce(rbind, Map(function(i) {list(clIds)}, seq_len(numOfMCMC - burnIn)))

centers <- vector(mode = "list", length = numOfTrees+1)
centers <- Reduce(rbind, Map(function(i) {list(centers)}, seq_len(numOfMCMC - burnIn)))

numOfBigNode <- rep(0, numOfMCMC - burnIn)
dd = matrix(0, nrow = numOfMCMC - burnIn, ncol = ncol(testData))
ssig = matrix(0, nrow = numOfMCMC - burnIn, ncol = ncol(testData))

subjAssignments <- matrix(0, nrow = numOfMCMC - burnIn, ncol = numOfSubjects)

## Conduct MCMC 
startTime <- Sys.time()

for (i in 0:numOfMCMC) {
  
  tssb$ResampleSubjectAssignments()

  tssb$ResampleSticksTrees()

  tssb$ResampleDataAssignments() 

  tssb$CullTree()

  tssb$ResampleNodeParameters()

  tssb$root$node$ResampleSigmaEach(Slice_sigma   = 1.0,
                                   Slice_stepOut = FALSE)

  tssb$root$node$ResampleHyperParams(Slice_sigma   = 0.1, 
                                     Slice_stepOut = FALSE)

  tssb$ResampleSticksZero(Slice_sigma   = 0.1,
                          Slice_stepOut = TRUE,
                          sd = NULL,
                          updateAlg = "slice")

  tssb$ResampleStickOrders(Slice_sigma   = 0.1,
                           Slice_stepOut = TRUE,
                           sd = NULL,
                           updateAlg = "slice",
                           Change = TRUE)

  tssb$SearchTreeYuan(Slice_sigma   = 0.1,
                      Slice_stepOut = TRUE,
                      sd = NULL,
                      updateAlg = "slice",
                      verbose = FALSE)

  tssb$ResampleBeta()

  tssb$KeepTrees()

  res_UnPost <- tssb$GetUnnormalizedPost()
  ll_UnPost[i] <- res_UnPost$ll
  
  tssb$oldSubjAssignments <- tssb$newSubjAssignments
  
  ## Save results after burn-in steps
  if (i > burnIn) { 
    ww <- tssb$GetMixture()
    subjAssignments[(i-burnIn), ] <- tssb$newSubjAssignments
    numOfBigNode[i-burnIn] <- sum(ww > 0.01)
    dd[i-burnIn, ] <- tssb$root$node$GetDrift()
    ssig[i-burnIn, ] <- tssb$root$node$GetSigma()
    
    ## Information for each subtree and tree0
    for (j in c(unique(tssb$newSubjAssignments), numOfTrees+1)) {
      if (j == (numOfTrees+1)) {
        resOfCl <- ClusterSummary_dphtssb_noinf(tssb$root, nrow(testData), tssbType = "yuan")
      } else {
        resOfCl <- ClusterSummary_dphtssb_noinf(tssb$root, nrow(testData), tssbType = "yuan", treeId = j)
      }
      eachCenter <- resOfCl$centers
      rownames(eachCenter) <- resOfCl$nodeIds
      
      clIds[[i-burnIn]][[j]] <- resOfCl$clIds
      centers[[i-burnIn]][[j]] <- eachCenter

      if (j == (numOfTrees+1)) {
        gg[[i-burnIn]][[j]] <- tssb$ConvertTssbToIgraph()$g
      } else {
        gg[[i-burnIn]][[j]] <- tssb$ConvertTssbToIgraph(j)$g
      }
    }
  }
  
  ## Output
  if (i <= burnIn) {
    if (i==0) cat(" Burn-in:") else if (i/10 == floor(i/10)) cat(paste0("+++", i))
  } else {
    if (i==burnIn+1) cat("\n MCMC sampling:") else if (i > burnIn & i/10 == floor(i/10)) cat(paste0("...", i))
  }
}

endTime <- Sys.time()
exeTime <- endTime - startTime


##########################
## Big Node Number Rule ##
##########################

## Select the iteration
uniqueBNN <- unique(numOfBigNode)
freq <- unlist(Map(function(i) {sum(numOfBigNode == uniqueBNN[i])}, 1:length(uniqueBNN)))
mostFreBNN <- uniqueBNN[which.max(freq)]
llAfterBurnIn <- ll_UnPost[(burnIn+1):numOfMCMC]
tssbNum <- which(llAfterBurnIn == max(llAfterBurnIn[numOfBigNode == mostFreBNN]))

## Final results
finalSubjAssignments <- subjAssignments[tssbNum, ]
finalClIds <- clIds[[tssbNum]][[numOfTrees + 1]]
finalSigma <- ssig[tssbNum, ]
finalDrift <- dd[tssbNum, ]

groups_subj <- sort(unique(finalSubjAssignments))
# sort() function is only for the convenience of the subsequent plotting step.

treeInd <- groups_subj[1]
ggFinal1 <- gg[[tssbNum]][[treeInd]]
subFinalClIds1 <- clIds[[tssbNum]][[treeInd]]
subFinalCenters1 <- centers[[tssbNum]][[treeInd]]

treeInd <- groups_subj[2]
ggFinal2 <- gg[[tssbNum]][[treeInd]]
subFinalClIds2 <- clIds[[tssbNum]][[treeInd]]
subFinalCenters2 <- centers[[tssbNum]][[treeInd]]

treeInd <- groups_subj[3]
ggFinal3 <- gg[[tssbNum]][[treeInd]]
subFinalClIds3 <- clIds[[tssbNum]][[treeInd]]
subFinalCenters3 <- centers[[tssbNum]][[treeInd]]


save(q0, tssb, testData, numOfSubjData, etaNormal, groups_subj, ggFinal1, ggFinal2, ggFinal3,
     subFinalClIds1, subFinalClIds2, subFinalClIds3,
     subFinalCenters1, subFinalCenters2, subFinalCenters3, 
     finalSigma, finalDrift, tssbNum, finalSubjAssignments, finalClIds, 
     file = "Real_result_iteration.RData")


