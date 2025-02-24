##############################################################################
## This code file is based on the original code file developed by
##
## Ke Yuan, Thomas Sakoparnig, Florian Markowetz and Niko Beerenwinkel. (2015) 
## BitPhylogeny: a probabilistic framework for reconstructing intra-tumor 
##    phylogenies. Genome Biology.
##
## 
## We changed the public member "sigma" to a private member and stored it
## in the root node. The priors of sigma and drift can be inverse Gamma 
## or uniform. We also added two private members "etaTheta" and "etaNormal", 
## which correspond to \eta_1 and \eta_2 in our model, respectively.
##############################################################################


#'@include node.R
NULL

#' R6 class for Normal node.
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords Normal node class
#' @field dataIds: data IDs
#' @field tssb: A TSSB object
#' @field children: descentant node objects
#' @field parent: parent node object
#' @method initialize
#' @method GetChildren
Normal <- R6::R6Class(
  classname = "Normal",
  inherit = Node,
  public = list(
    params = NULL,

    initialize = function(parent = NULL, tssb = NULL, dataDims = 1, depth = 0, numOfTrees = NULL,
                          drift = 1, 
                          etaNormal = 1,
                          etaTheta  = 1,
                          priorDriftMin = 0.01,
                          priorDriftMax = 1,
                          priorDriftScale = 1,
                          priorDriftShape = 1,
                          priorSigmaMin = 0.01,
                          priorSigmaMax = 1,
                          priorSigmaScale = 1,
                          priorSigmaShape = 1,
                          priorEtaNormalA = 0.1,
                          priorEtaNormalB = 1,    # eta_normal ~ Unif(A, B)
                          priorEtaThetaA  = 0.1,
                          priorEtaThetaB  = 1,    # eta_theta  ~ Unif(A, B)
                          initMean = array(0, dim=c(1, dataDims)),
                          priorTypeSigma = c("invGamma", "unif"),
                          priorTypeDrift = c("invGamma", "unif")) {
 
      super$initialize(parent = parent, tssb = tssb, dataDims = dataDims, depth = depth, numOfTrees = numOfTrees)
      if (is.null(parent)) {
        private$drift <- if (priorTypeDrift == "invGamma") {
          rinvgamma(dataDims, shape = priorDriftShape, scale = priorDriftScale)
          } else if (priorTypeDrift == "unif") {
          runif(dataDims, min = priorDriftMin, max = priorDriftMax)
          }  # row vec
        private$sigma <- if (priorTypeSigma == "invGamma") {
          rinvgamma(dataDims, shape = priorSigmaShape, scale = priorSigmaScale)
          } else if (priorTypeSigma == "unif") {
          runif(dataDims, min = priorSigmaMin, max = priorSigmaMax)
          }  # row vec
        private$etaNormal <- etaNormal
        private$etaTheta  <- etaTheta
        private$priorDriftMin   <- priorDriftMin
        private$priorDriftMax   <- priorDriftMax
        private$priorDriftScale <- priorDriftScale
        private$priorDriftShape <- priorDriftShape
        private$priorSigmaMin   <- priorSigmaMin
        private$priorSigmaMax   <- priorSigmaMax
        private$priorSigmaScale <- priorSigmaScale
        private$priorSigmaShape <- priorSigmaShape
        private$priorEtaNormalA <- priorEtaNormalA
        private$priorEtaNormalB <- priorEtaNormalB
        private$priorEtaThetaA  <- priorEtaThetaA
        private$priorEtaThetaB  <- priorEtaThetaB
        private$initMean <- initMean
        self$params <- rnorm(dataDims, mean=initMean, sd = sqrt(private$drift))
        private$priorTypeSigma <- priorTypeSigma
        private$priorTypeDrift <- priorTypeDrift
      } else {
        drift <- self$GetDrift()
        etaTheta  <- self$GetEtaTheta()
        self$params <- rnorm(dataDims, mean=parent$params, sd = sqrt(etaTheta^(self$depth) * drift))
      }
    },

    # Drift
    GetPriorDriftShape = function() {
      if (is.null(private$parent)) {
        private$priorDriftShape
      } else {
        private$parent$GetPriorDriftShape()
      }
    },
    
    GetPriorDriftScale = function() {
      if (is.null(private$parent)) {
        private$priorDriftScale
      } else {
        private$parent$GetPriorDriftScale()
      }
    },
    
    GetPriorDriftMin = function() {
      if (is.null(private$parent)) {
        private$priorDriftMin
      } else {
        private$parent$GetPriorDriftMin()
      }
    },
    
    GetPriorDriftMax = function() {
      if (is.null(private$parent)) {
        private$priorDriftMax
      } else {
        private$parent$GetPriorDriftMax()
      }
    },
    
    GetDrift = function() {
      if (is.null(private$parent)) {
        private$drift
      } else {
        private$parent$GetDrift()
      }
    },
    
    # sigma
    GetPriorSigmaShape = function() {
      if (is.null(private$parent)) {
        private$priorSigmaShape
      } else {
        private$parent$GetPriorSigmaShape()
      }
    },
    
    GetPriorSigmaScale = function() {
      if (is.null(private$parent)) {
        private$priorSigmaScale
      } else {
        private$parent$GetPriorSigmaScale()
      }
    },
    
    GetPriorSigmaMin = function() {
      if (is.null(private$parent)) {
        private$priorSigmaMin
      } else {
        private$parent$GetPriorSigmaMin()
      }
    },
    
    GetPriorSigmaMax = function() {
      if (is.null(private$parent)) {
        private$priorSigmaMax
      } else {
        private$parent$GetPriorSigmaMax()
      }
    },
    
    GetSigma = function() {
      if (is.null(private$parent)) {
        private$sigma
      } else {
        private$parent$GetSigma()
      }
    },
    
    GetPriorTypeSigma = function() {
      if (is.null(private$parent)) {
        private$priorTypeSigma
      } else {
        private$parent$GetPriorTypeSigma()
      }
    },
    
    # eta
    GetEtaNormal = function() {
      if (is.null(private$parent)) {
        private$etaNormal
      } else {
        private$parent$GetEtaNormal()
      }
    },
    
    GetEtaTheta = function() {
      if (is.null(private$parent)) {
        private$etaTheta
      } else {
        private$parent$GetEtaTheta()
      }
    },

    GetLogProb = function(x) {
      etaNormal <- self$GetEtaNormal()
      sigma <- self$GetSigma()
      depth <- self$depth
      if (!is.null(dim(x))) { 
        sum(apply(x, MARGIN = 1, function(datum) {dnorm(datum, mean = self$params, sd = sqrt(etaNormal^(depth) * sigma), log = T)} ))
      } else {
        sum(dnorm(x, mean = self$params, sd=sqrt(etaNormal^(depth) * sigma), log = T))
      }
    },

    GetNodeLogProb = function(treeId) {
      if (!missing(treeId)) {
        self$GetLogProb(self$GetData(treeId))
      } else {
        self$GetLogProb(self$GetData())
      }
    },
    
    GetParamLogPrior = function() {
      etaTheta <- self$GetEtaTheta()
      mean <- if (is.null(private$parent)) private$initMean else private$parent$params
      res <- sum(dnorm(self$params, mean = mean, sd = sqrt(etaTheta^(self$depth) * self$GetDrift()), log = T))
      return(res)
    },
    
    
    ## Update 
    ResampleParams = function() {
      nodeData <- self$GetData()  
      drift <- self$GetDrift()
      sigma <- self$GetSigma()
      etaNormal <- self$GetEtaNormal()
      etaTheta  <- self$GetEtaTheta()
      depth <- self$depth
      
      numOfData <- self$GetNumOfLocalData() 
      numOfChildren <- length(private$children)  ## num of children: W_ch

      if (numOfData == 0) {
        dataMean = 0
      } else if (numOfData == 1){
        dataMean = nodeData
      } else {
        dataMean = colMeans(nodeData)
      }

      if (is.null(private$parent)) {
        parentParams <- private$initMean
      } else {
        parentParams <- private$parent$params
      }

      if (numOfChildren ==0 ) {  ## num of children is always W
        childParamsMean = 0
      } else {
        childParams <- Reduce(
          rbind,
          Map(function(x) {x$params}, private$children),
          c()
        )
        childParamsMean <- colMeans(childParams)
      }

      # Construct prior for node mean
      priorParamsMean <- (etaTheta*parentParams + numOfChildren*childParamsMean)/(numOfChildren + etaTheta)
      priorParamsCov <- etaTheta^(depth+1)*drift / (numOfChildren + etaTheta)

      # Posterior for node mean
      if (numOfData == 0) {
        postParamsMean <- priorParamsMean
        postParamsCov <- priorParamsCov
      } else {
        postParamsCov <- (priorParamsCov^(-1) + numOfData/(etaNormal^(depth) * sigma))^(-1)  
        postParamsMean <- (priorParamsMean/priorParamsCov +
                             numOfData*dataMean/(etaNormal^(depth) * sigma)) * postParamsCov
      }
      self$params <- rnorm(n = self$dataDims,
                           mean = postParamsMean,
                           sd = sqrt(postParamsCov))

      invisible(self)
    },
    
    
    # Update sigma for each
    ResampleSigmaEach = function(Slice_sigma = 1.0, Slice_stepOut = T) {
      etaNormal <- private$etaNormal
      priorSigmaScale <- private$priorSigmaScale
      priorSigmaShape <- private$priorSigmaShape
      
      if(!is.null(private$parent)) {
        stop("Can only update sigma from root!")
      }
      
      ComputeDataLlh <- function(eachSigma, num) {
        
        if (private$priorTypeSigma == "invGamma") {
          if (eachSigma <= 0) { 
            return (-Inf)
          }
        }
        
        if (private$priorTypeSigma == "unif") {
          if (eachSigma < self$GetPriorSigmaMin() || eachSigma > self$GetPriorSigmaMax()) {
            return (-Inf)
          }
        }
        
        
        Descend <- function(root, depth = 0) {
          numOfLocalData <- root$GetNumOfLocalData()
          localData <- root$GetData()                
          if (numOfLocalData == 0) {
            llh = 0
          } else if (numOfLocalData == 1){
            llh = dnorm(localData[num], mean = root$params[num], sd=sqrt(etaNormal^(depth) * eachSigma), log = TRUE)
          } else {
            llh = sum(dnorm(localData[, num], mean = root$params[num], sd=sqrt(etaNormal^(depth) * eachSigma), log = TRUE))
          }
          
          children <- root$GetChildren()
          for (i in seq_along(children)) {
            child <- children[[i]]
            llh <- llh + Descend(child, depth = depth+1)
          }
          return(llh)
        }
        
        res <- Descend(self)
        if (private$priorTypeSigma == "invGamma") {
          res <- res + log(dinvgamma(eachSigma, shape = priorSigmaShape, scale = priorSigmaScale))
        }
        return(res)
      }
      private$sigma <- unlist(Map(function(j) {SliceSampler(private$sigma[j], ComputeDataLlh, sigma = Slice_sigma, stepOut = Slice_stepOut, num = j)}, 1:self$dataDims)) 
    },
    
    
    
    ResampleSigmaEach_v2 = function() {
      if(!is.null(private$parent)) {
        stop("Can only update sigma from root!")
      }
      
      etaNormal <- private$etaNormal
      priorSigmaScale <- private$priorSigmaScale
      priorSigmaShape <- private$priorSigmaShape
      
      roots <- self$tssb$GetRootList()
      numOfAllData <- sum(self$tssb$numOfSubjData)
      postSigmaShape <- priorSigmaShape + numOfAllData / 2
      
      tmpSum <- Reduce(
        sum,
        Map(function(i) {
          node <- roots[[i]]$node
          if (node$GetNumOfLocalData()) {
            localData <- node$GetData()
            if (!is.null(dim(localData))) {  
              rowSums(apply(localData, 1, function(datum) {(datum - node$params)^2})) / etaNormal^(node$depth)
            } else {
              (localData - node$params)^2 / etaNormal^(node$depth)
            }
            
          } else {
            0
          }
        },
        seq_along(roots)
        ),
        0
      )
      
      postSigmaScale <- priorSigmaScale + tmpSum / 2
      private$sigma <- rinvgamma(self$dataDims, shape = postSigmaShape, scale = postSigmaScale)
      
    },
    
    
    

    # Update drift
    ResampleHyperParams = function(Slice_sigma = 1.0, Slice_stepOut = T) {
      etaTheta <- self$GetEtaTheta()

      if(!is.null(private$parent)) {
        stop("Can only update hypers from root!")
      }
      
      ComputeDriftLlh <- function(eachDrift, num) {
        if (private$priorTypeDrift == "invGamma") {
          if (eachDrift <= 0) { 
            return (-Inf)
          }
        }
        
        if (private$priorTypeDrift == "unif") {
          if (eachDrift < self$GetPriorDriftMin() || eachDrift > self$GetPriorDriftMax()) {
            return (-Inf)
          }
        }
        
        Descend <- function(root, depth = 1) {
          llh <- 0
          children <- root$GetChildren()
          for (i in seq_along(children)) {
            child <- children[[i]]
            llh <- llh + dnorm(child$params[num], root$params[num], sqrt(etaTheta^(depth)*eachDrift), log = TRUE)
            llh <- llh + Descend(child, depth = depth+1)
          }
          return(llh)
        }
        
        res <- Descend(self) + dnorm(self$params[num], private$initMean[num], sqrt(eachDrift), log = TRUE)
        if (private$priorTypeDrift == "invGamma") {
          res <- res + log(dinvgamma(eachDrift, private$priorDriftShape, private$priorDriftScale)) 
        }
        return(res)
      }
      private$drift <- unlist(Map(function(j) {SliceSampler(private$drift[j], ComputeDriftLlh, sigma = Slice_sigma, stepOut = Slice_stepOut, num = j)}, 1:self$dataDims)) 
    },
    
    
    # Update two Eta
    ResampleHyperEta = function(sampleEtaNormal = T,
                                sampleEtaTheta  = T,
                                Slice_sigma = 1.0, 
                                Slice_stepOut = T) {
      if(!is.null(private$parent)) {
        stop("Can only update hypers from root!")
      }
      priorEtaNormalA <- private$priorEtaNormalA
      priorEtaNormalB <- private$priorEtaNormalB
      priorEtaThetaA  <- private$priorEtaThetaA
      priorEtaThetaB  <- private$priorEtaThetaB
      drift <- private$drift
      sigma <- private$sigma
      
      # Compute log-likelihood for "etaNormal"
      ComputeEtaNormalLlh <- function(etaNormal) {
        if (etaNormal < priorEtaNormalA || etaNormal > priorEtaNormalB) {
          return(-Inf)
        }
        Descend <- function(root, depth = 0) {
          llh <- 0
          children <- root$GetChildren()
          for (i in seq_along(children)) {
            child <- children[[i]]
            numOfLocalData <- child$GetNumOfLocalData()  
            localData <- child$GetData()  
            
            if (numOfLocalData == 0) {
              llh = llh + 0
            } else if (numOfLocalData == 1){
              llh = llh + sum(dnorm(localData, mean = child$params, sd=sqrt(etaNormal^(depth+1) * sigma), log = TRUE))
            } else {
              llh = llh + sum(apply(localData, MARGIN = 1, function(datum) {
                dnorm(datum, mean = child$params, sd = sqrt(etaNormal^(depth+1) * sigma), log = TRUE)
              } ))
            }
            
            llh <- llh + Descend(child, depth = depth+1)
          }
          return(llh)
        }
        return(Descend(self))
      }
      
      # Compute log-likelihood for "etaTheta"
      ComputeEtaThetaLlh <- function(etaTheta) {
        if (etaTheta < priorEtaThetaA || etaTheta > priorEtaThetaB) {
          return(-Inf)
        }
        Descend <- function(root, depth = 0) {
          llh <- 0
          children <- root$GetChildren()
          for (i in seq_along(children)) {
            child <- children[[i]]
            llh <- llh + sum(dnorm(child$params, root$params, sqrt(etaTheta^(depth+1) * drift), log = TRUE))
            llh <- llh + Descend(child, depth = depth+1)
          }
          return(llh)
        }
        return(Descend(self))
      }
      
      if (sampleEtaNormal) {
        private$etaNormal <- SliceSampler(private$etaNormal, ComputeEtaNormalLlh, sigma = Slice_sigma, stepOut = Slice_stepOut)
      }
      if (sampleEtaTheta) {
        private$etaTheta  <- SliceSampler(private$etaTheta,  ComputeEtaThetaLlh,  sigma = Slice_sigma, stepOut = Slice_stepOut)
      }
      
    }
    
  ),
  private = list(
    priorTypeDrift = NULL,
    priorTypeSigma = NULL,
    drift = NA,
    sigma = NA,
    etaNormal = NA,
    etaTheta = NA,
    priorDriftMin = NA,
    priorDriftMax = NA,
    priorDriftScale = NA,
    priorDriftShape = NA,
    priorSigmaMin = NA,
    priorSigmaMax = NA,
    priorSigmaScale = NA,
    priorSigmaShape = NA,
    priorEtaNormalA = NA,
    priorEtaNormalB = NA,
    priorEtaThetaA  = NA,
    priorEtaThetaB  = NA,
    initMean = NA
  )
)