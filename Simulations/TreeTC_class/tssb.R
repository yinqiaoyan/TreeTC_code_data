##############################################################################
## This code file is based on the original code file developed by
##
## Ke Yuan, Thomas Sakoparnig, Florian Markowetz and Niko Beerenwinkel. (2015) 
## BitPhylogeny: a probabilistic framework for reconstructing intra-tumor 
##    phylogenies. Genome Biology.
##
## 
## We added public member "csum" to avoid redundant computation. 
## We added public methods: 
##    * GetUnnormalizedPost(): Compute the unnormalized posterior proportional 
##      to the first equation in Section "Posterior sampling".
##    * GetMixtureFromPath(): Compute the weight by path and treeId directly.
##    * FindNodeWeights(): Determine all the elements in the finite set
##      mentioned in sampling procedure (d).
##    * GetRootList(): Get the root list in tree0.
##    * GetPathsAll(): Get the paths of nodes in tree0.
##    * ResetPaths(): Reorganize paths after some updates.
##############################################################################


#'TSSB is a R6 object of tree-structured stick breaking
#' process
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords data
#' @field data
#' @field dpAlpha
#' @field dpGamma
#' @field dpLambda
#' @field minDepth
#' @field maxDepth
#' @field root
#' @field assignments
#' @field numOfData
#' @field minDpAlpha
#' @field maxDpAlpha
#' @field minDpLambda
#' @field minDpLambda
#' @field minDpGamma
#' @field maxDpGamma
#' @method new TSSB constuctor
#' @method FindNode Find node in tree, generate new node as needed
#' @method GetMixture Compute the mixing weights and get the corresponding nodes
#' @method ConvertTssbToIgraph Convert the tssb root list to an igraph object
#' @method CullTree Remove empty leaf nodes
#' @method GetLogMarginalDataLikelihood Compute the log marginal data likelihood
TSSB <- R6::R6Class(
  classname = "TSSB",
  public = list(
    data = NULL,                
    root = NULL,         
    assignments = NULL,         
    oldSubjAssignments = NULL,  # new added! numeric vec, len = numOfSubjects
    newSubjAssignments = NULL,  # new added! numeric vec, len = numOfSubjects
    numOfSubjects = NA,         # new added! double
    numOfSubjData = NULL,       # new added! numeric vec, len = numOfSubjects
    numOfTrees = NA,            # new added! double
    beta = NULL,                # new added! numeric vec, len = numOfTrees
    dpKappa = NA,               # new added! double
    # Fields ------------------------------------------------------------------
    dpAlpha0 = NA,       
    dpGamma0 = NA,
    dpLambda0 = NA,
    dpAlpha1 = NA,              # new added!
    dpGamma1 = NA,              # new added!
    minDepth = NA,      
    maxDepth = NA,
    minDpAlpha0 = NA,
    maxDpAlpha0 = NA,
    minDpLambda0 = NA,
    maxDpLambda0 = NA,
    minDpGamma0 = NA,
    maxDpGamma0 = NA,
    minDpAlpha1 = NA,           # new added!
    maxDpAlpha1 = NA,           # new added!
    minDpGamma1 = NA,           # new added!
    maxDpGamma1 = NA,           # new added!
    csum = NA,                  # new added!
    
    # Methods -----------------------------------------------------------------
    initialize = function(rootNode = emptyenv(),
                          dpKappa = 1,           # new added!
                          dpAlpha0 = 1,
                          dpGamma0 = 1,
                          dpLambda0 = 1,
                          dpAlpha1 = 1,          # new added!
                          dpGamma1 = 1,          # new added!
                          data = NULL,           # Input：data matrix!
                          numOfSubjects = NULL,  # new added!
                          numOfSubjData = NULL,  # new added!
                          minDepth = 0,
                          maxDepth = 25,
                          minDpAlpha0 = 1e-6, maxDpAlpha0 = 50,
                          minDpLambda0 = 1e-6, maxDpLambda0 = 1,
                          minDpGamma0 = 1e-6, maxDpGamma0 = 50,
                          minDpAlpha1 = 1e-6, maxDpAlpha1 = 50,  # new added!
                          minDpGamma1 = 1e-6, maxDpGamma1 = 50,  # new added!
                          initType = c("unif", "beta")
                          ) {
      if (!is.environment(rootNode) ||
            identical(rootNode, emptyenv()) ||
            !"Node" %in% class(rootNode)) {
        stop("Root node must be specified")  
      }

      # init some numbers
      self$numOfSubjects <- numOfSubjects
      self$numOfSubjData <- numOfSubjData
      self$csum <- c(0, cumsum(self$numOfSubjData))
      self$numOfTrees <- rootNode$numOfTrees
      self$newSubjAssignments <- rep(0, self$numOfSubjects)
      self$data <- data          
      self$dpKappa  <- dpKappa            # new added!
      self$dpAlpha0 <- dpAlpha0
      self$dpGamma0 <- dpGamma0
      self$dpAlpha1 <- dpAlpha1           # new added!
      self$dpGamma1 <- dpGamma1           # new added!
      self$minDepth <- minDepth
      self$maxDepth <- maxDepth
      self$dpLambda0 <- dpLambda0
      self$minDpAlpha0 <- minDpAlpha0
      self$maxDpAlpha0 <- maxDpAlpha0
      self$minDpLambda0 <- minDpLambda0
      self$maxDpLambda0 <- maxDpLambda0
      self$minDpGamma0 <- minDpGamma0
      self$maxDpGamma0 <- maxDpGamma0
      self$minDpAlpha1 <- minDpAlpha1    # new added!
      self$maxDpAlpha1 <- maxDpAlpha1    # new added!
      self$minDpGamma1 <- minDpGamma1    # new added!
      self$maxDpGamma1 <- maxDpGamma1    # new added!
      self$root <- list(
        node        = rootNode,                                                    
        keep0       = NULL,     # (bool vec, len = numOfChildren)                       Record if children are active in tree0
        keepTrees   = NULL,     # (bool mat, nrow = numOfTrees, ncol = numOfChildren)   Record if children are active in each subtree
        main0       = NULL,     # (double)                                              Record main in tree0
        mainTrees   = NULL,     # (double vec, len = numOfTrees)                        Record main in each subtree
        sticks0     = NULL,     # (double vec, len = numOfChildren)                     Record sticks in tree0
        sticksTrees = NULL,     # (double mat, nrow = numOfTrees, ncol = numOfChildren) Record sticks in each subtree
        children    = list())
      newMain0 <- if (self$minDepth == 0) BoundBeta(1, 1, self$dpAlpha0) else 0
      self$root$main0 <- newMain0
      self$root$mainTrees <- if (self$minDepth == 0) {
        BoundBeta(self$numOfTrees, dpAlpha1 * newMain0, dpAlpha1 * (1 - newMain0))
        } else {rep(0, self$numOfTrees)}
      rootNode$tssb <- self 
      rootNode$path <- c()   
      
      
      # init subject assignments -----------------------------------------------
      if (initType == "beta") {
        DPsticks <- c(BoundBeta(self$numOfTrees - 1, 1, self$dpKappa), 1)
        edges <- SticksToEdges(DPsticks)
        self$beta <- diff(c(0, edges))
        self$oldSubjAssignments <- sample(1:self$numOfTrees,   
                                          self$numOfSubjects,   
                                          replace = T,
                                          prob = self$beta)
      } else if (initType == "unif") {
        self$beta <- rep(1/self$numOfTrees, self$numOfTrees)
        self$oldSubjAssignments <- sample(1:self$numOfTrees,    
                                          self$numOfSubjects,   
                                          replace = T,
                                          prob = self$beta)
      }
      
      
      
      # Draw data assignments and a tree in init step! -------------------------
      k <- 0
      for (j in 1:self$numOfSubjects) {
        treeId <- self$oldSubjAssignments[j]
        numOfData <- self$numOfSubjData[j]
        uSamples <- runif(numOfData)
        for (n in (k+1):(k+numOfData)) { 
          res <- self$FindNode(uSamples[n-k], treeId)  
          self$assignments <- c(self$assignments, res$node)          
          res$node$AddDatum(n)           
          res$node$AddDatum(n, treeId)   
          self$root <- res$root 
        }
        k <- k + numOfData
      }
      
    },
    
    
    

    FindNode = function(u, treeId) {
      
      Descend <- function(root, u, path = c(), depth = 0, flag = 0) {
        if (u < root$mainTrees[treeId]) {
          return(list(node = root$node, path = path, root = root, flag = flag))
        } else {
          u <- (u - root$mainTrees[treeId]) / (1 - root$mainTrees[treeId])
          
          while (length(root$children) == 0  
                 || u > (1 - prod(1 - root$sticksTrees[treeId,]))  
                 ) {
            flag <- 1
            newStick0 <- BoundBeta(1, 1, self$dpGamma0)
            root$sticks0 <- c(root$sticks0, newStick0) 
            newStickTrees <- BoundBeta(self$numOfTrees, self$dpGamma1 * newStick0, self$dpGamma1 * (1-newStick0))
            root$sticksTrees <- cbind(root$sticksTrees, as.matrix(newStickTrees)) 
            newChild <- list(
              list(
                node = root$node$Spawn(), 
                keep0 = NULL,
                keepTrees = NULL, 
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
          }
          
          edges <- c(0, SticksToEdges(root$sticksTree[treeId,]))
          index <- sum(u > edges)
          u     <- (u - edges[index]) / (edges[index+1] - edges[index])
          res <- Descend(root$children[[index]], u, path, depth+1, flag)  
          node <- res$node  
          path <- c(index, res$path)  
          root$children[[index]] <- res$root  

          flag = res$flag
          return(list(node = node, path = path, root = root, flag = flag))
          }
      }
      return(Descend(self$root, u))
    },
    

    
    FindNodeWeights = function(u_ij, treeId) {
      
      Descend <- function(root, u_ij, mass, depth = 0) {
        weight <- mass * root$mainTrees[treeId] 
        node <- if (weight >= u_ij) root$node else NULL
        rest <- mass * (1.0 - root$mainTrees[treeId])
        
        if (rest < u_ij) {
          return(list(node = node, root = root))
        } else {
          edges <- SticksToEdges(root$sticksTrees[treeId,])
          while (length(root$children) == 0  
                 || rest * (1 - edges[length(edges)]) >= u_ij  
          ) {
            newStick0 <- BoundBeta(1, 1, self$dpGamma0)
            root$sticks0 <- c(root$sticks0, newStick0) 
            newStickTrees <- BoundBeta(self$numOfTrees, self$dpGamma1 * newStick0, self$dpGamma1 * (1-newStick0))
            root$sticksTrees <- cbind(root$sticksTrees, as.matrix(newStickTrees)) 
            newChild <- list(
              list(
                node = root$node$Spawn(), 
                keep0 = NULL,
                keepTrees = NULL, 
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
            
            edges <- SticksToEdges(root$sticksTrees[treeId,])
          }
          
          weights <- diff(c(0, edges))
          for (i in seq_along(root$children)) {
            res <- Descend(root$children[[i]], u_ij, mass*(1.0-root$mainTrees[treeId])*weights[i], depth+1)
            node <- c(node, res$node)
            root$children[[i]] <- res$root
          }
          
          return(list(node = node, root = root))
          
        }
      }
      return(Descend(self$root, u_ij, 1.0))
    },
    
    

    GetNumOfDataAllNodes = function() {
      Descend <- function(root) {
        numOfData <- root$node$GetNumOfLocalData()
        for (i in seq_along(root$children)) {
          numOfData <- c(numOfData, Descend(root$children[[i]]))
        }
        return(numOfData)
      }
      return(Descend(self$root))
    },
    
    

    # Get the weights of all nodes
    GetMixture = function(treeId) {
      if (!missing(treeId)) {
        Descend <- function(root, mass) {
          weight <- mass * root$mainTrees[treeId] 
          edges <- SticksToEdges(root$sticksTrees[treeId,])
          weights <- diff(c(0, edges))
          
          for (i in seq_along(root$children)) {
            child <- root$children[[i]]
            weight <- c(weight, Descend(child, mass*(1.0-root$mainTrees[treeId])*weights[i]))
          }
          return(weight = weight)
        }
      } else {
        Descend <- function(root, mass) {
          weight <- mass * root$main0 
          edges <- SticksToEdges(root$sticks0)
          weights <- diff(c(0, edges))
          
          for (i in seq_along(root$children)) {
            child <- root$children[[i]]
            weight <- c(weight, Descend(child, mass*(1.0-root$main0)*weights[i]))
          }
          return(weight = weight)
        }
      }
      
      return(Descend(self$root, 1.0))
    },
    
    
    
    GetMixtureFromPath = function(path, treeId) {
      Descend <- function(root, mass, depth = 0) {
        weight <- mass * root$mainTrees[treeId] 
        edges <- SticksToEdges(root$sticksTrees[treeId,])
        weights <- diff(c(0, edges))
        
        if (depth == length(path)) return(weight = weight)
        
        i <- path[depth+1]
        child <- root$children[[i]]
        weight <- Descend(child, mass*(1.0-root$mainTrees[treeId])*weights[i], depth + 1)
        return(weight = weight)
      }
      
      return(Descend(tssb$root, 1.0))
    },
    
    
    # Get the paths of nodes in tree0
    GetPathsAll = function() {
      Descend <- function(root) {
        pathVec <- list(root$node$path)
        for (i in seq_along(root$children)) {
          pathVec <- rbind(pathVec, Descend(root$children[[i]]))
        }
        return(pathVec)
      }
      return(Descend(self$root))
    },
    

    ResetPaths = function() {
      Descend <- function(root) {
        for (i in seq_along(root$children)) {
          child <- root$children[[i]]
          child$node$path <- c(root$node$path, i)
          root$children[[i]] <- Descend(child)
        }
        return(root)
      }
      self$root <- Descend(self$root)
    },
    
    
    
    # Get root list 
    GetRootList = function() {
      Descend <- function(root) {
        rootList <- list(root) 
        for (i in seq_along(root$children)) {
          rootList <- c(rootList, Descend(root$children[[i]]))
        }
        return(rootList)
      }
      return(Descend(self$root))
    },
    
    
    
    
    GetNumOfActiveNodes = function() {
      Descend <- function(root) {
        num = 0
        keep <- root$keep0
        if (sum(keep) < 1) {
          return(num)
        }
        num <- num + sum(keep)
        for (i in which(keep)) {
          num <- num + Descend(root$children[[i]])
        }
        return(num)
      }
      return(Descend(self$root) + 1)  # +1 root node
    },
    
    
   

    ConvertTssbToIgraph = function(treeId) { 
      if (missing(treeId)) {
        edges <- SticksToEdges(self$root$sticks0)
        weights <- diff(c(0, edges))
        g <- igraph::graph.empty(directed = TRUE)
        g <- g + igraph::vertex(name = "X",
                                size = self$root$node$GetNumOfLocalData())
        
        Descend <- function(root, name, mass, g) {
          if (length(root$sticks0) < 1){
            return(list(total = mass, g=g))
          } else {
            total <- 0
            edges <- SticksToEdges(root$sticks0)
            weights <- diff(c(0, edges))
            
            for (i in 1:length(root$children)) {
              child <- root$children[[i]]
              childName <- paste(name, i, sep = "-")
              childMass <- mass * weights[i] * child$main0
              g <- g + igraph::vertex(name = childName,
                                      size = child$node$GetNumOfLocalData())
              g <- g + igraph::edge(name, childName,
                                    Value = child$node$GetNumOfLocalData())
              tmp <- Descend(child,
                             childName,
                             mass*weights[i]*(1.0 - child$main0),
                             g)
              g <- tmp$g
              total = total + childMass + tmp$total
            }
            return(list(total=total, g=g))
          }
        }
      } else {
        edges <- SticksToEdges(self$root$sticksTrees[treeId,])
        weights <- diff(c(0, edges))
        g <- igraph::graph.empty(directed = TRUE)
        g <- g + igraph::vertex(name = "X",
                                size = self$root$node$GetNumOfLocalData(treeId))
        
        Descend <- function(root, name, mass, g) {
          if (length(root$sticksTrees[treeId,]) < 1){
            return(list(total = mass, g=g))
          } else {
            total <- 0
            edges <- SticksToEdges(root$sticksTrees[treeId,])
            weights <- diff(c(0, edges))
            
            for (i in 1:length(root$children)) {
              child <- root$children[[i]]
              childName <- paste(name, i, sep = "-")
              childMass <- mass * weights[i] * child$mainTrees[treeId]
              g <- g + igraph::vertex(name = childName,
                                      size = child$node$GetNumOfLocalData(treeId))
              g <- g + igraph::edge(name, childName,
                                    Value = child$node$GetNumOfLocalData(treeId))
              tmp <- Descend(child,
                             childName,
                             mass*weights[i]*(1.0 - child$mainTrees[treeId]),
                             g)
              g <- tmp$g
              total = total + childMass + tmp$total
            }
            return(list(total=total, g=g))
          }
        }
      }
      res = Descend(self$root, "X", 1, g)
      return(res)
    },
    
    
    
    ConvertSubtreeToIgraph = function(treeId) {  
      if (missing(treeId)) {
        edges <- SticksToEdges(self$root$sticks0)
        weights <- diff(c(0, edges))
        g <- igraph::graph.empty(directed = TRUE)
        g <- g + igraph::vertex(name = "X",
                                size = self$root$node$GetNumOfLocalData())
        
        Descend <- function(root, name, mass, g) {
          keepChildren <- root$keep0
          if (sum(keepChildren) < 1){  
            return(list(total = mass, g=g))
          } else {
            total <- 0
            edges <- SticksToEdges(root$sticks0[keepChildren])
            weights <- diff(c(0, edges))
            
            for (i in which(keepChildren)) {
              child <- root$children[[i]]
              childName <- paste(name, i, sep = "-")
              childMass <- mass * weights[i] * child$main0
              g <- g + igraph::vertex(name = childName,
                                      size = child$node$GetNumOfLocalData())
              g <- g + igraph::edge(name, childName,
                                    Value = child$node$GetNumOfLocalData())
              tmp <- Descend(child,
                             childName,
                             mass*weights[i]*(1.0 - child$main0),
                             g)
              g <- tmp$g
              total = total + childMass + tmp$total
            }
            return(list(total=total, g=g))
          }
        }
      } else {
        edges <- SticksToEdges(self$root$sticksTrees[treeId,])
        weights <- diff(c(0, edges))
        g <- igraph::graph.empty(directed = TRUE)
        g <- g + igraph::vertex(name = "X",
                                size = self$root$node$GetNumOfLocalData(treeId))
        
        Descend <- function(root, name, mass, g) {
          keepChildren <- root$keepTrees[treeId,]
          if (sum(keepChildren) < 1){  
            return(list(total = mass, g=g))
          } else {
            total <- 0
            edges <- SticksToEdges(root$sticksTrees[treeId,][keepChildren])
            weights <- diff(c(0, edges))
            
            for (i in which(keepChildren)) {
              child <- root$children[[i]]
              childName <- paste(name, i, sep = "-")
              childMass <- mass * weights[i] * child$mainTrees[treeId]
              g <- g + igraph::vertex(name = childName,
                                      size = child$node$GetNumOfLocalData(treeId))
              g <- g + igraph::edge(name, childName,
                                    Value = child$node$GetNumOfLocalData(treeId))
              tmp <- Descend(child,
                             childName,
                             mass*weights[i]*(1.0 - child$mainTrees[treeId]),
                             g)
              g <- tmp$g
              total = total + childMass + tmp$total
            }
            return(list(total=total, g=g))
          }
        }
      }
      res = Descend(self$root, "X", 1, g)
      return(res)
    },
    
    
    
    CullTree = function() {
      Descend <- function(root) {
        res <- unlist(Map(Descend, root$children), recursive = F, use.names = F)

        if (length(res) == 0) {  
          return(list(counts = root$node$GetNumOfLocalData(), root = root))
        }
        counts <- unlist(res[seq(1, length(res), by = 2)])
        root$children <- res[seq(2, length(res), by = 2)]  
        keep <- which(counts != 0)

        sapply(root$children[which(counts == 0)], function(x) x$node$Kill())

        if (length(keep) == 0) {
          root$sticks0 <- NULL
          root$sticksTrees <- NULL
          root$children <- NULL
        } else {
          root$sticks0 <- root$sticks0[keep]
          root$sticksTrees <- as.matrix(root$sticksTrees[, keep])  
          root$children <- root$children[keep]
        }
        return(list(counts = sum(counts) + root$node$GetNumOfLocalData(),
                    root = root))
      }
      res <- Descend(self$root)
      self$root <- res$root
      self$ResetPaths()  
      invisible(self)
    },
    
    
    
    KeepTrees = function() {
      Descend <- function(root) {
        res <- unlist(Map(Descend, root$children), recursive = F, use.names = F)
        
        if (length(res) == 0) {  
          root$keep0 <- NULL
          root$keepTrees <- NULL
          return(list(counts0 = root$node$GetNumOfLocalData(),
                      countsTrees = unlist(Map(root$node$GetNumOfLocalData, 1:self$numOfTrees)),  
                      root = root))
        }
        counts0            <- unlist(res[seq(1, length(res), by = 3)])  
        countsTrees <- Reduce(cbind, res[seq(2, length(res), by = 3)], c())  
        root$children             <- res[seq(3, length(res), by = 3)]   
        
        root$keep0 <- counts0 != 0
        root$keepTrees <- countsTrees != 0
        
        return(list(counts0 = sum(counts0) + root$node$GetNumOfLocalData(),
                    countsTrees = rowSums(countsTrees) + unlist(Map(root$node$GetNumOfLocalData, 1:self$numOfTrees)),
                    root = root))  
      }
      self$root <- Descend(self$root)$root
      invisible(self)
    },

    
    
   
    GetLogMarginalDataLikelihood = function() {
      ll.output <- 0
      roots <- self$GetRootList()
      for (treeId in unique(self$newSubjAssignments)) {
        weights <- self$GetMixture(treeId)
        
        ll <- Reduce(
          sum,
          Map(function(i) {
            node <- roots[[i]]$node
            if (node$GetNumOfLocalData(treeId)) {
              node$GetNumOfLocalData(treeId)*weights[i] + node$GetNodeLogProb(treeId)
            } else {
              0
            }
          },
          seq_along(weights)
          ),
          0
        )
        
        ll.output <- ll.output + ll
      }
      return(list(ll = ll.output, nn = roots)) 
    },
    
    
    
    # Compute the unnormalized posterior of (C, S, theta)
    # Proportion to the first equation in Section "Posterior sampling"
    GetUnnormalizedPost = function() {
      ll.output <- 0
      roots <- self$GetRootList()
      seq_roots <- seq_along(roots)
      for (treeId in unique(self$newSubjAssignments)) {
  
        weights <- self$GetMixture(treeId)
        beta_treeId <- self$beta[treeId]
        
        ll <- Reduce(
          sum,
          Map(function(i) {
            node <- roots[[i]]$node
            if (node$GetNumOfLocalData(treeId)) {
              return(node$GetNumOfLocalData(treeId)*log(weights[i]) + 
                       node$GetNodeLogProb(treeId) +
                       node$GetParamLogPrior()) 
            } else { 
              return(0)
            }
          },
          seq_roots
          ),
          0
        )
        ll.output <- ll.output + ll + 
          sum(self$newSubjAssignments == treeId) * log(beta_treeId) 
      }
      
      
      return(list(ll = ll.output, nn = roots))
    }
  )
)



