##############################################################################
## This code file is based on the original code file developed by
##
## Ke Yuan, Thomas Sakoparnig, Florian Markowetz and Niko Beerenwinkel. (2015) 
## BitPhylogeny: a probabilistic framework for reconstructing intra-tumor 
##    phylogenies. Genome Biology.
## 
##
## We added two public members: "numOfTrees" and "path", and 
## replaced "dataIds" with "dataIdsAll" and "dataIdsTrees" 
## to let this code file be used for multiple tree structures in our model.
## The public methods:
## HasData(), AddDatum(), RemoveDatum(), GetNumOfLocalData(), 
## GetNumOfSubTreeData()
## are also replaced with multiple tree structure versions.
##############################################################################


#' Node is a R6 object of each cluster in TSSB
#'
#'
#' @docType class
#' @importFrom R6 R6Class
#' @export
#' @format An \code{\link{R6Class}} generator object
#' @keywords Node class
#' @field dataIds: data IDs
#' @field tssb: A TSSB object
#' @field children: descentant node objects
#' @field parent: parent node object
#' @method new
#' @method GetChildren
#' @method GetParent
#' @method SetParent
#' @method AddChild
#' @method RemoveChild
#' @method Kill
#' @method Spawn lay an egg
#' @method HasData
#' @method AddDatum
#' @method RemoveDatum
#' @method GetNumOfLocalData
#' @method GetNumOfSubTreeData
#' @method GetData
#' @method GetAncestors

Node <- R6::R6Class(
  classname = "Node",

  public = list(
    # Fields ------------------------------------------------------------------
    dataIdsAll = NULL,
    dataIdsTrees = NULL,
    tssb = NULL,
    dataDims = NULL,
    depth = NULL,
    numOfTrees = NULL,
    path = NULL,        

    # Methods -----------------------------------------------------------------
    initialize = function(parent = NULL, 
                          tssb = NULL, 
                          dataDims = NULL,
                          depth = NULL,
                          numOfTrees = NULL) {
      if (!is.null(parent)) {
        private$parent <- parent
        parent$AddChild(self)
      } else {
        private$parent <- parent
      }
      if (!is.null(tssb) ) self$tssb <- tssb
      
      self$numOfTrees <- numOfTrees
      self$dataIdsTrees <- vector("list", length = self$numOfTrees)
      self$dataDims <- dataDims
      self$depth <- depth
      },


    GetChildren = function() {
      private$children
    },

    GetParent = function() {
      private$parent
    },

    SetParent = function(parent = NULL) {
      private$parent <- parent
    },

    AddChild = function(child) {
      if (!missing(child)) {
        child$SetParent(self)
        private$children <- c(private$children, child)
      }
      invisible(self)
    },

    RemoveChild = function(child) {
      if (!missing(child)) {
        private$children <- unlist(Filter(
          Negate(
            function(x) identical(child, x)
            ),
          private$children))
      }
      invisible(self)
    },

    Kill = function() {
      if (!is.null(private$parent))  {
        private$parent$RemoveChild(self)
      }
      private$children = NULL
      private$parent = NULL
    },

    Spawn = function() {
      return(get(class(self)[1])$new(parent = self, 
                                     tssb = self$tssb, 
                                     dataDims = self$dataDims,
                                     depth = self$depth+1,
                                     numOfTrees = self$numOfTrees)) 
    },

    HasData = function(treeId) {
      if (!missing(treeId)) {
        if ( length(self$dataIdsTrees[[treeId]])>0 ) {
          return(TRUE)
        } else {
          return(
            sum(
              unlist(sapply(private$children, function(x) x$HasData(treeId)))
            )
            > 0)
        }
      } else {
        if ( length(self$dataIdsAll)>0 ) {
          return(TRUE)
        } else {
          return(
            sum(
              unlist(sapply(private$children, function(x) x$HasData()))
            )
            > 0)
        }
      }
      
    },

    AddDatum = function(id, treeId) {
      if (!missing(id)) { 
        if (!missing(treeId)) {
          if (!any(self$dataIdsTrees[[treeId]]==id)) self$dataIdsTrees[[treeId]] <- c(self$dataIdsTrees[[treeId]], id)
        } else {
          if (!any(self$dataIdsAll==id)) self$dataIdsAll <- c(self$dataIdsAll, id)
        }
      }
    },

    RemoveDatum = function(id, treeId) {
      if (!missing(id)) {  
        if (!missing(treeId)) {
          if (!id %in% self$dataIdsTrees[[treeId]]) {
            warning("id is not found in dataIds, nothing is removed")
          } else {
            self$dataIdsTrees[[treeId]] <- self$dataIdsTrees[[treeId]][self$dataIdsTrees[[treeId]] != id]
          }
        } else {
          if (!id %in% self$dataIdsAll) {
            warning("id is not found in dataIds, nothing is removed")
          } else {
            self$dataIdsAll <- self$dataIdsAll[self$dataIdsAll != id]
          }
        }
      }
    },

    GetNumOfLocalData = function(treeId) {
      if (!missing(treeId)) length(self$dataIdsTrees[[treeId]]) else length(self$dataIdsAll)
    },

    GetNumOfSubTreeData = function(treeId) {
      if (!missing(treeId)) {
        Reduce(
          sum,
          Map(function(x) x$GetNumOfSubTreeData(treeId),
              private$children),
          length(self$dataIdsTrees[[treeId]])
        )
      } else {
        Reduce(
          sum,
          Map(function(x) x$GetNumOfSubTreeData(),
              private$children),
          length(self$dataIdsAll)
        )
      }
    },

    GetData = function(treeId) {
      if (!missing(treeId)) {
        if (dim(self$tssb$data)[2] > 1) {
          return(self$tssb$data[self$dataIdsTrees[[treeId]],])
        } else {
          return(as.matrix(self$tssb$data[self$dataIdsTrees[[treeId]],]))
        }
      } else {
        if (dim(self$tssb$data)[2] > 1) {
          return(self$tssb$data[self$dataIdsAll,])
        } else {
          return(as.matrix(self$tssb$data[self$dataIdsAll,]))
        }
      }
    },

    GetAncestors = function() {
      ancestors = c()
      if (is.null(private$parent)) {
        return(list(self))
      } else {
        ancestors <- c(private$parent$GetAncestors(), self)
        return(ancestors)
      }
    }

    ), # end of public

  private = list(
    # Fields
    children = NULL,
    parent = NULL
    )
  )




