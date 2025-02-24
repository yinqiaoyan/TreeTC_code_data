##############################################################################
## This code file is based on the original code file developed by
##
## Ke Yuan, Thomas Sakoparnig, Florian Markowetz and Niko Beerenwinkel. (2015) 
## BitPhylogeny: a probabilistic framework for reconstructing intra-tumor 
##    phylogenies. Genome Biology.
##
## 
## We added some public methods: 
##    * ResampleSubjectAssignments(): Update the source group membership S^{(j)}
##    * ResampleDataAssignments(): Update the observation group membership 
##      C_{i}^{(j)}. To prevent an infinite for-loop due to underflow, 
##      we set a lower bound for u_{ij}.
##    * ResampleSticksTrees(): Update sticks \nu^{(l)} and \psi^{(l)}
##    * ResampleSticksZero(): Update sticks \nu and \psi
##    * ResampleBeta(): Update the allocation weights \omega_{l}
##############################################################################


#'@include tssb.R
NULL

#' R6 class for inference via MCMC for TSSB.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data
#' @method ResampleAssignments Sample data assignments with Adams's slice-retrospective sampler
#' @method ResampleSticks Sample sticks
#' @method ResampleStickOrders Sample stick orders
#' @method ResampleNodeParameters Sample node parameters
#' @method ResampleHypers Sample TSSB hyper parameters

TssbMCMC <- R6::R6Class(
  classname = "TssbMCMC",
  inherit = TSSB,
  public = list(
    
    initialize = function(...) {
      super$initialize(...)  
    },
    
    
    ResampleSubjectAssignments = function() {
      csum <- self$csum
      weightsMat <- Reduce(rbind, Map(self$GetMixture, 1:self$numOfTrees))
      roots <- self$GetRootList()
      pathsAll <- self$GetPathsAll()  # Return：list of paths of each root
      numOfCurrentRoots <- length(roots)
      seq_numOfCurrentRoots <- seq_len(numOfCurrentRoots)  

      for (j in 1:self$numOfSubjects) {
        sumOfLogProb <- 0
        sumOfLogWeights <- rep(0, self$numOfTrees)
        pathsOfData <- Map(function(i) {self$assignments[[i]]$path}, (csum[j] + 1):csum[j + 1])  # Return: list
        uniqPaths <- unique(pathsOfData)  # Return: list of unique elements
        for (n in seq_along(uniqPaths)) {
          dataIds <- Filter(function(i) identical(uniqPaths[[n]], pathsOfData[[i]]), seq_len(csum[j+1]-csum[j])) + csum[j]
          nodeId  <- Filter(function(i) identical(uniqPaths[[n]], pathsAll[[i]]), seq_numOfCurrentRoots)
          node <- roots[[nodeId]]$node
          
          sumOfLogProb <- sumOfLogProb + node$GetLogProb(self$data[dataIds,])
          sumOfLogWeights <- sumOfLogWeights + length(dataIds) * log(weightsMat[, nodeId])
        }
        postWeights <- log(self$beta) + sumOfLogWeights + sumOfLogProb
        postWeights <- postWeights - max(postWeights)
        postWeights <- exp(postWeights) / sum(exp(postWeights))  # log-sum-exp trick
        self$newSubjAssignments[j] <- sample(1:self$numOfTrees, 1, prob = postWeights)
      }
    },
    
    
    
    ResampleDataAssignments = function() {
      csum <- self$csum
      for (j in 1:self$numOfSubjects) {
        treeIdNew <- self$newSubjAssignments[j]
        treeIdOld <- self$oldSubjAssignments[j]

        for (n in (csum[j]+1):csum[j+1]) {
          
          indices = self$assignments[[n]]$path
          weightOld <- self$GetMixtureFromPath(indices, treeIdNew)
          u_ij <- runif(1, min = 0, max = weightOld)
          if (u_ij < 1e-18) u_ij = 1e-18
          
          res <- self$FindNodeWeights(u_ij, treeIdNew)
          self$root <- res$root
          nodeList <- res$node
          
          if (is.list(nodeList)) {
            postWeights <- unlist(Map(function(node) {node$GetLogProb(self$data[n,])}, nodeList))
            postWeights <- postWeights - max(postWeights)
            postWeights <- exp(postWeights) / sum(exp(postWeights))  # log-sum-exp trick
            newNode <- nodeList[[sample(seq_along(nodeList), 1, prob = postWeights)]]
          } else {
            newNode <- nodeList
          }
          
          
          if (!identical(newNode, self$assignments[[n]]) || treeIdNew != treeIdOld) {  
            self$assignments[[n]]$RemoveDatum(n)
            self$assignments[[n]]$RemoveDatum(n, treeIdOld)
            newNode$AddDatum(n)
            newNode$AddDatum(n, treeIdNew)
            self$assignments[[n]] <- newNode
          }
        }
      }
      invisible(self)
    },
    
    
    
    # resample sticks nu^l, psi^l for pi^l
    ResampleSticksTrees = function() {
      Descend <- function(root, treeId, depth=0) {
        dataDown <- 0
        indices <- seq_along(root$children)  
        
        for (i in sort(indices, decreasing = T)) {
          res <- Descend(root$children[[i]], treeId, depth+1)
          childData <- res$nodeData
          root$children[[i]] <- res$root
          postAlpha <- self$dpGamma1 * root$sticks0[i] + childData
          postBeta <- self$dpGamma1 * (1 - root$sticks0[i]) + dataDown
          root$sticksTrees[treeId, i] <- BoundBeta(1, postAlpha, postBeta)
          dataDown <- dataDown + childData
        }
        
        dataHere <- root$node$GetNumOfLocalData(treeId) 
        postAlpha <- self$dpAlpha1 * root$main0 + dataHere
        postBeta <- self$dpAlpha1 * (1 - root$main0) + dataDown
        root$mainTrees[treeId] <- if (self$minDepth > depth) {0.0
        } else if (self$maxDepth > depth) {
          BoundBeta(1, postAlpha, postBeta)
        } else {1.0}
        return(list(nodeData = dataHere + dataDown, root = root))
      }
      
      for (l in 1:self$numOfTrees) {
        self$root <- Descend(self$root, l)$root
      }
      invisible(self)
    },
    
    
    
    # resample sticks nu, psi for pi^0
    ResampleSticksZero = function(Slice_sigma = 1.0, 
                                  Slice_stepOut = T,
                                  sd = 0.01, 
                                  updateAlg = c("slice", "mh")) {
      ComputeMain0Llh <- function(main0, root) {
        if (main0 < 0 || main0 > 1) {
          return(-Inf)
        }
        llh <- sum(dbeta(root$mainTrees, self$dpAlpha1*main0, self$dpAlpha1*(1-main0), log = TRUE))
        llh <- llh + dbeta(main0, 1, self$dpAlpha0 * self$dpLambda0^(root$node$depth), log = TRUE)
        return(llh)
      }

      ComputeStick0Llh <- function(stick0, root, num) {
        if (stick0 < 0 || stick0 > 1) {
          return(-Inf)
        }
        llh <- sum(dbeta(root$sticksTrees[, num], self$dpGamma1*stick0, self$dpGamma1*(1-stick0), log = TRUE))
        llh <- llh + dbeta(stick0, 1, self$dpGamma0, log = TRUE)
        return(llh)
      }
      
      if (updateAlg == "mh") {
        rProposal <- function(initX, sd) {
          return(rnorm(1, mean = initX, sd = sd))
        }
        
        dLogProposal <- function(propX, initX, sd) {
          return(dnorm(propX, mean = initX, sd = sd, log = TRUE))
        }
      }

      Descend <- function(root) {
        root$main0 <- if (updateAlg == "slice") {
          SliceSampler(root$main0,
                       ComputeMain0Llh,
                       sigma = Slice_sigma,
                       stepOut = Slice_stepOut,
                       root = root)
        } else if (updateAlg == "mh") {
          MH(root$main0,
             ComputeMain0Llh,
             dLogProposal,
             rProposal,
             sd = sd,
             root = root)
        }
          
        for (i in seq_along(root$children)) {
          root$sticks0[i] <- if (updateAlg == "slice") {
            SliceSampler(root$sticks0[i],
                         ComputeStick0Llh,
                         sigma = Slice_sigma,
                         stepOut = Slice_stepOut,
                         root = root,
                         num  = i)
          } else if (updateAlg == "mh") {
            MH(root$sticks0[i],
               ComputeStick0Llh,
               dLogProposal,
               rProposal,
               sd = sd,
               root = root,
               num  = i)
          }
            
          if (length(root$node$path) < self$maxDepth - 1) {
            root$children[[i]] <- Descend(root$children[[i]])
          }
        }
        return(root = root)
      }

      self$root <- Descend(self$root)
    },
    
    
    # Resample stickOrder in tree0
    ResampleStickOrders = function(Slice_sigma = 1.0, 
                                   Slice_stepOut = T,
                                   sd = 0.01, 
                                   updateAlg = c("slice", "mh"),
                                   Change = T) {
      Descend <- function(root, depth = 0) {
        if (length(root$children)==0) {
          return(root)
        }
        newOrder <- c()
        represented <- Filter(function(i) root$children[[i]]$node$HasData(),
                              seq_along(root$children)) 
        allWeights <- diff(c(0,SticksToEdges(root$sticks0)))

        while (length(represented)!=0) {
          u <- runif(1)
          while (TRUE) {
            subIndices <- Filter(function(x) ! x %in% newOrder, seq_along(root$sticks0))
            subWeights <- c(allWeights[subIndices], 1-sum(allWeights))
            subWeights <- subWeights/sum(subWeights)
            index <- sum(u > cumsum(subWeights))
            if (index == length(subIndices)) {
              newStick0 <- BoundBeta(1, 1, self$dpGamma0)
              root$sticks0 <- c(root$sticks0, newStick0) 
              newStickTrees <- BoundBeta(self$numOfTrees, self$dpGamma1 * newStick0, self$dpGamma1 * (1-newStick0))
              root$sticksTrees <- cbind(root$sticksTrees, as.matrix(newStickTrees))  
              newChild <- list(
                list(
                  node = root$node$Spawn(), 
                  keep = NULL,  
                  main0 = NULL,
                  mainTrees = NULL,
                  sticks0 = NULL,
                  sticksTrees = NULL,
                  children = NULL
                )
              )
              
              newChild[[1]]$node$path <- c(root$node$path, length(root$children)+1)
              newMain0 <- if (self$minDepth > (depth+1)) {0.0
              } else if (self$maxDepth > (depth+1)) {
                BoundBeta(1, 1, (self$dpLambda0^(depth+1))*self$dpAlpha0)
              } else {1.0}
              newChild[[1]]$main0 <- newMain0
              newChild[[1]]$mainTrees <- if (self$minDepth > (depth+1)) {
                rep(0.0, self$numOfTrees)  
              } else if (self$maxDepth > (depth+1)) {
                BoundBeta(self$numOfTrees, self$dpAlpha1 * newMain0, self$dpAlpha1 * (1 - newMain0))
              } else {rep(1.0, self$numOfTrees)}
              
              root$children <- c(root$children, newChild)  
              allWeights <- diff(c(0, SticksToEdges(root$sticks0))) 
            } else {
              index <- subIndices[index+1]
              break
            }
          }
          newOrder <- c(newOrder, index)
          represented <- represented[represented != index]
        }

        newChildren <- c()
        if (Change) {
          for (k in newOrder) {
            root$children[[k]] <- Descend(root$children[[k]], depth + 1)
            child <- root$children[k] 
            newChildren <- c(newChildren, child)  
          }
        } else {
          for (k in newOrder) {
            child <- root$children[k]
            newChildren <- c(newChildren, child) 
            root$children[[k]] <- Descend(root$children[[k]], depth + 1)
          }
        }
        

        lapply(Filter(function(x) ! x %in% newOrder, seq_along(root$sticks0)),
               function(k) {root$children[[k]]$node$Kill()
                            root$children[[k]] <- list()})

        root$children <- newChildren
        root$sticks0 <- BoundBeta(length(newChildren), 1, self$dpGamma0)
        root$sticksTrees <- matrix(0, nrow = self$numOfTrees, ncol = length(newChildren))

        return(root)
      }
      self$root <- Descend(self$root)
      self$ResampleSticksTrees()
      self$ResampleSticksZero(Slice_sigma = Slice_sigma, 
                              Slice_stepOut = Slice_stepOut,
                              sd = sd,
                              updateAlg = updateAlg)  
      self$ResetPaths()  
      invisible(self)
    },

    
    # Resample node params theta and sigma
    ResampleNodeParameters = function() {
      Descend <- function(root) {
        Map(Descend, root$children)
        root$node$ResampleParams()
        return(root=root)
      }
      self$root <- Descend(self$root)  
      invisible(self)
    },
    
    

    ResampleHypers = function(sampleDpAlpha0 = T,
                              sampleDpLambda0 = T,
                              sampleDpGamma0 = T,
                              sampleDpAlpha1 = T,
                              sampleDpGamma1 = T) {

      # Compute all the log-likelihoods ----------------------------------------
      # Compute llh for main0
      ComputeDpMain0Llh <- function(params,
                                    sample = c("alpha0", "lambda0"),
                                    fixParams = NULL) {
        if (sample == "alpha0") {
          dpAlpha0 <- params
          dpLambda0 <- fixParams
        } else if (sample == "lambda0") {
          dpAlpha0 <- fixParams
          dpLambda0 <- params
        } else {
          stop("sample has to be either 'alpha' or 'lambda' ")
        }

        if (dpAlpha0 < self$minDpAlpha0 || dpAlpha0 > self$maxDpAlpha0) {
          return(-Inf)
        }
        if (dpLambda0 < self$minDpLambda0 || dpLambda0 > self$maxDpLambda0) {
          return(-Inf)
        }
        Descend <- function(root, depth = 0) {
          llh <- if (self$minDepth <= depth && self$maxDepth > depth) { 
            dbeta(root$main0, 1, dpLambda0^depth*dpAlpha0, log = TRUE)
          } else {0}

          for (i in seq_along(root$children)) {
            res <- Descend(root$children[[i]], depth + 1)
            llh <- llh + res$llh
            root$children[[i]] = res$root
          }
          return(list(llh = llh, root = root))
        }
        return(Descend(self$root)$llh)
      }

      # Compute llh for sticks0
      ComputeDpBranch0Llh <- function(dpGamma0) {
        if (dpGamma0 < self$minDpGamma0 || dpGamma0 > self$maxDpGamma0) {
          return(-Inf)
        }
        Descend <- function(root) {
          llh <- 0
          for (i in seq_along(root$children)) {
            llh <- llh + dbeta(root$sticks0[i], 1, dpGamma0, log = TRUE)
            res <- Descend(root$children[[i]])
            llh <- llh + res$llh
            root$children[[i]] <- res$root
          }
          return(list(llh = llh, root = root))
        }
        return(Descend(self$root)$llh)
      }
      
      # Compute llh for mainTrees
      ComputeDpMainTreesLlh <- function(dpAlpha1) {
        if (dpAlpha1 < self$minDpAlpha1 || dpAlpha1 > self$maxDpAlpha1) {
          return(-Inf)
        }
        Descend <- function(root, depth = 0) {
          llh <- if (self$minDepth <= depth && self$maxDepth > depth) { 
            sum(dbeta(root$mainTrees, dpAlpha1 * root$main0, dpAlpha1 * (1-root$main0), log = TRUE))
          } else {0}
          for (i in seq_along(root$children)) {
            res <- Descend(root$children[[i]], depth + 1)
            llh <- llh + res$llh
            root$children[[i]] = res$root
          }
          return(list(llh = llh, root = root))
        }
        return(Descend(self$root)$llh)
      }
      
      # Compute llh for sticksTrees
      ComputeDpBranchTreesLlh <- function(dpGamma1) {
        if (dpGamma1 < self$minDpGamma1 || dpGamma1 > self$maxDpGamma1) {
          return(-Inf)
        }
        Descend <- function(root) {
          llh <- 0
          for (i in seq_along(root$children)) {
            llh <- llh + sum(dbeta(root$sticksTrees[, i], dpGamma1*root$sticks0[i], dpGamma1*(1-root$sticks0[i]), log = TRUE))
            res <- Descend(root$children[[i]])
            llh <- llh + res$llh
            root$children[[i]] <- res$root
          }
          return(list(llh = llh, root = root))
        }
        return(Descend(self$root)$llh)
      }
      
      # Apply slice sampler ----------------------------------------------------

      if (sampleDpAlpha0) {
        self$dpAlpha0 <- SliceSampler(self$dpAlpha0,
                                     ComputeDpMain0Llh,
                                     sample = "alpha0",
                                     fixParams = self$dpLambda0)
      }
      if (sampleDpLambda0) {
        self$dpLambda0 <- SliceSampler(self$dpLambda0,
                                      ComputeDpMain0Llh,
                                      sample = "lambda0",
                                      fixParams = self$dpAlpha0)

      }
      if (sampleDpGamma0) {
        self$dpGamma0 <- SliceSampler(self$dpGamma0, ComputeDpBranch0Llh)
      }
      if (sampleDpAlpha1) {
        self$dpAlpha1 <- SliceSampler(self$dpAlpha1, ComputeDpMainTreesLlh)
      }
      if (sampleDpGamma1) {
        self$dpGamma1 <- SliceSampler(self$dpGamma1, ComputeDpBranchTreesLlh)
      }
      invisible(self)
    },
    
    
    # Resample outer DP weights beta
    ResampleBeta = function() {
      subjNums <- unlist(Map(function(i) {sum(self$newSubjAssignments == i)}, 1:self$numOfTrees))
      postDPsticks <- unlist(Map(function(l) {
        BoundBeta(1, 1 + subjNums[l], self$dpKappa + sum(subjNums[(l+1):self$numOfTrees]))
      }, 1:(self$numOfTrees-1)))
      postDPsticks <- c(postDPsticks, 1)
      edges <- SticksToEdges(postDPsticks)
      self$beta <- diff(c(0, edges))
    },
    
    
    # swap-nodes step
    SearchTreeYuan = function(Slice_sigma = 1.0, 
                              Slice_stepOut = T,
                              sd = 0.01, 
                              updateAlg = c("slice", "mh"), 
                              verbose = F) {
      res <- self$GetUnnormalizedPost()
      post <- rexp(1) + res$ll
      roots <- res$nn

      len <- length(roots)
      
      if (len > 1) {
        nodeAnum <- sample(1:len, 1)
        nodeBnum <- sample(1:len, 1)
        while (nodeAnum == nodeBnum) nodeBnum <- sample(1:len, 1)
        
        SwapNodes <- function(nodeAnum, nodeBnum) {
     
          # swap nodeParams and dataIds
          nodeA <- roots[[nodeAnum]]$node
          nodeB <- roots[[nodeBnum]]$node
          
          paramsA <- nodeA$params
          dataIdsAllA <- nodeA$dataIdsAll
          dataIdsTreesA <- nodeA$dataIdsTrees
          
          nodeA$params <- nodeB$params
          Map(function(i) {nodeA$RemoveDatum(i)}, nodeA$dataIdsAll)
          Map(function(treeId) {
            Map(function(i) {nodeA$RemoveDatum(i, treeId)}, nodeA$dataIdsTrees[[treeId]])
            }, 1:self$numOfTrees)
          Map(function(i) {
            nodeA$AddDatum(i)
            self$assignments[[i]] <- nodeA
            }, nodeB$dataIdsAll)
          Map(function(treeId) {
            Map(function(i) {nodeA$AddDatum(i, treeId)}, nodeB$dataIdsTrees[[treeId]])
          }, 1:self$numOfTrees)
          
          nodeB$params <- paramsA
          
          Map(function(i) {nodeB$RemoveDatum(i)}, nodeB$dataIdsAll)
          Map(function(treeId) {
            Map(function(i) {nodeB$RemoveDatum(i, treeId)}, nodeB$dataIdsTrees[[treeId]])
          }, 1:self$numOfTrees)
          Map(function(i) {
            nodeB$AddDatum(i)
            self$assignments[[i]] <- nodeB
            }, dataIdsAllA)
          Map(function(treeId) {
            Map(function(i) {nodeB$AddDatum(i, treeId)}, dataIdsTreesA[[treeId]])
          }, 1:self$numOfTrees)
          
          # swap mains
          main0A <- roots[[nodeAnum]]$main0
          mainTreesA <- roots[[nodeAnum]]$mainTrees
          main0B <- roots[[nodeBnum]]$main0
          mainTreesB <- roots[[nodeBnum]]$mainTrees
          
          Descend <- function(root, nodeNum) {
            if (nodeNum == nodeAnum) {
              root$main0 <- main0B
              root$mainTrees <- mainTreesB
            }
            if (nodeNum == nodeBnum) {
              root$main0 <- main0A
              root$mainTrees <- mainTreesA
            }
            for (i in seq_along(root$children)) {
              nodeNum <- nodeNum + 1
              res <- Descend(root$children[[i]], nodeNum)
              root$children[[i]] <- res$root
              nodeNum <- res$nodeNum
            }
            return(list(root=root, nodeNum=nodeNum))
          }
          
          if (main0A != 1 && main0B != 1) { 
            self$root <- Descend(root=self$root, nodeNum=1)$root
          }
        }
        
        SwapNodes(nodeAnum, nodeBnum)
        self$ResampleSticksTrees()
        self$ResampleSticksZero(Slice_sigma = Slice_sigma, 
                                Slice_stepOut = Slice_stepOut,
                                sd = sd,
                                updateAlg = updateAlg)
        
        post_new <- self$GetUnnormalizedPost()$ll
        if (post_new < post) {
          SwapNodes(nodeAnum, nodeBnum)
          self$ResampleSticksTrees()
          self$ResampleSticksZero(Slice_sigma = Slice_sigma, 
                                  Slice_stepOut = Slice_stepOut,
                                  sd = sd,
                                  updateAlg = updateAlg)
        } else {if (verbose) cat("Successful swap!\n")}
      }
    }
  )
  
)





