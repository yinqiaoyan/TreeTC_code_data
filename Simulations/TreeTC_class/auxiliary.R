##############################################################
################## Some auxiliary functions ##################
##############################################################


# MCMC Summary for DP-HTSSB
ClusterSummary_dphtssb_noinf = function(tssbRoot, numOfAllData, tssbType, treeId) {
  dims <- length(tssbRoot$node$params)
  
  if (missing(treeId)) {  # For tree0
    Descend <- function(root, clIds=array(0, dim = numOfAllData), 
                        centers=NULL, nodeIds=c(),
                        name, nodeNum=1) {
      dataIds <- root$node$dataIdsAll
      params  <- root$node$params
      
      nodeIds <- c(nodeIds, name)
      if (length(dataIds) > 0) {
        clIds[dataIds] <- nodeNum
        centers <- rbind(centers, c(params, 1))  # 1 means this center has data
        nodeNum <- nodeNum + 1
      } else {
        centers <- rbind(centers, c(params, 0))  # 0 means no data
      }
      if (tssbType == "tssbdw") children <- root$children[root$keep]  
      if (tssbType == "yuan") children <- root$children
      for (i in seq_along(children)) {
        childName <- paste(name, i, sep = "-")
        res <- Descend(children[[i]], clIds, 
                       centers, nodeIds, 
                       childName, nodeNum)
        clIds <- res$clIds
        centers <- res$centers
        nodeIds <- res$nodeIds
        nodeNum <- res$nodeNum
      }
      return(list(clIds = clIds, 
                  centers = centers,
                  nodeIds = nodeIds,
                  nodeNum = nodeNum))
    }
    
    res = Descend(tssbRoot, name="X")
    return(res)
  } else {  # For subtree treeId
    Descend <- function(root, clIds=array(0, dim = numOfAllData), 
                        centers=NULL, nodeIds=c(),
                        name, nodeNum=1) {
      dataIds <- root$node$dataIdsTrees[[treeId]] 
      params  <- root$node$params
      
      nodeIds <- c(nodeIds, name)
      if (length(dataIds) > 0) {
        clIds[dataIds] <- nodeNum
        centers <- rbind(centers, c(params, 1))
        nodeNum <- nodeNum + 1
      } else {
        centers <- rbind(centers, c(params, 0))
      }
      if (tssbType == "tssbdw") children <- root$children[root$keep]  
      if (tssbType == "yuan") children <- root$children
      for (i in seq_along(children)) {
        childName <- paste(name, i, sep = "-")
        res <- Descend(children[[i]], clIds, 
                       centers, nodeIds, 
                       childName, nodeNum)
        clIds <- res$clIds
        centers <- res$centers
        nodeIds <- res$nodeIds
        nodeNum <- res$nodeNum
      }
      return(list(clIds = clIds, 
                  centers = centers,
                  nodeIds = nodeIds,
                  nodeNum = nodeNum))
    }
    
    res = Descend(tssbRoot, name="X")
    res$clIds = res$clIds[res$clIds != 0]
    return(res)
  }
}