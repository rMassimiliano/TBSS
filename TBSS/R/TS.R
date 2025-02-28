#'General TS class for TBSS
#'@slot data contains the data in a \code{data.frame} one row for each leaf of the tree.  The number of columns varies depending on the analysis, while the first column is the name of a leaf. 
#'@slot tree a\code{data.frame} with two columns representing the hierarchical structure in a edge list format. Each row represents a directional edge between a node and its parent node, so the two columns are a nodeID and parentID. The IDs are treated as \code{character}.
#'@slot leaves is a vector with the names of the leaves, typically the first column of the slot \code{data}
#'@slot nodes a vector with the names of all the nodes in the tree, i.e. the first column of \code{tree}
#'@slot mapNodesLeaves is a named list. Each element of the list is named after a node, while the content of the list indicates the leaves associated with that node. 
#'@slot nodesToTest is a vector defining the nodes to be tested. Some nodes could be excluded from the testing procedure  (e.g., because there are too few data points) and can be excluded from this list.
#'@slot nodeSS is a \code{data.frame} reporting the sufficient statistics for each node in the tree. The number of columns varies depending on the analysis, and the first column is the name of a node. 
#'@slot LRT is a named vector containing the node-specific discrepancy measure computed on the observed data for each node in \code{nodesToTest}.
#' @slot B is the number of MC replicates used to approximate the null distribution of LRT
#'@slot LRT_H0 is a \code{matrix} with the MC-distribution of LRT. Each column is an MC replicate.
#' @name TS
#' @rdname TS
#' @exportClass TS
setClass(Class = "TS", slots = c(data = 'data.frame',
                                 tree ='data.frame',
                                 leaves = 'character',
                                 nodes = 'character',
                                 mapNodesLeaves = 'list',
				 nodesToTest ='character',
                                 nodeSS = 'data.frame',
                                 LRT ='numeric',
                                 B = 'numeric',
                                 LRT_H0 = 'matrix'))


#' Populate the leaves slot
#' @description Populate the slot leaves of TS with a vector of the unique values of the first column of data (leaves ID). 
#' @param object an object of class TS
listLeaves = function(object)
{
 object@leaves =  unique(object@data[,1])
 return(object)
}

#' Populate the slot nodes
#' @description Populate the slot node of TS with a vector containing the unique value of the first column of the slot \code{tree}. Empty \code{character} are not allowed
#' @param object an object of class TS
listNodes = function(object)
{
  object@nodes =  setdiff(unique(unlist(object@tree)),"")
  names(object@nodes) = NULL
  object
}


#' Populate the slot mapNodesLeaves
#' @description Create a named list. Each element of the list is named after a node and contains the ID of leaves linked to that node. The function returns an object of class TS with a populated slot \code{mapNodeLeaves}.
#' @param object an object of class TS
#' @param prune if \code{prune = TRUE} remove nodes for which we do not have data. If nodes are removed, print a warning. 
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.

mapNodesLeaves = function(object, prune = TRUE,parallel = FALSE, ncpus = NULL)
{
  if(length(object@leaves) == 0) object = listLeaves(object)
  if(length(object@nodes) == 0) object = listNodes(object)
  get_path = function(node, prune)
  {
     ## Leaves are only connected to
     if(node %in% object@leaves)
     {
       return(node)
     }
     tmp = getLeaf(node,object@tree,object@leaves)
     tmp = tmp[which(tmp %in% object@leaves)]
     if(prune)
     {
      if(length(tmp) > 0){return(tmp)}
      else{return(NULL)}
     }
     else
     {
       return(tmp)
     }
  }
  if(parallel)
  {
   suppressMessages(require(parallel))
   suppressMessages(require(TBSS))
   cl <- makeCluster(ncpus)
    clusterExport(cl, "get_path", envir = environment())
    clusterExport(cl, "getLeaf", envir = environment())
    clusterExport(cl, "object", envir = environment())
    clusterExport(cl, "prune", envir = environment())

    nodeList = parLapply(cl, object@nodes,\(x) get_path(x,prune))
    ind = parLapplyLB(cl,nodeList,length)
    ind = which(ind > 0)
    node_removed = length(nodeList) - length(ind)
    stopCluster(cl)
   }
   else
   {
   nodeList = lapply(object@nodes, \(x) get_path(x,prune))
   ind = which(sapply(nodeList,length) > 0)
   node_removed = length(nodeList) - length(ind)
   }

   names(nodeList)= object@nodes
   nodeList = nodeList[ind]
  object@mapNodesLeaves = nodeList
  if(node_removed>0) warning(sprintf("We removed %i nodes from the analysis because no data have been provided for the corresponding leaves", node_removed))
  return(object)
}


#' Compute node-specific sufficient statistics
#' @description generic function to compute node-specific sufficient statistics. The function populates the slot \code{nodeSS}. 
#' @param object an object of class TS
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.

setGeneric("computeSS", function(object, parallel = FALSE, ncpus = NULL) {
  standardGeneric("computeSS")
})

setMethod("computeSS", signature(object = "TS"), function(object,parallel = FALSE, ncpus = NULL) {
  stop("To compute node specific sufficient statistics we need a model. Use or implement an appropriate function.")
})




#' Compute node-specific discrepancy measures
#' @description generic function to compute node-specific discrepancy measures. The function populates the slot \code{nodeSS}. 
#' @param object an object of class TS
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.

setGeneric("computeLRT", function(object, parallel = FALSE, ncpus = NULL) {
  standardGeneric("computeLRT")
})



# Compute the LRT. The expression depends on the model
setMethod("computeLRT", signature(object = "TS"), function(object,parallel = FALSE, ncpus = NULL) {
  stop("The LRT formula depends on the type of analysis. Use or implement an appropriate class.")
})

#' Generate data under the null distribution
#' @description generic function to generate data under the null distribution. Instances of this function are used in \code{monteCarlo} 
#' @param object an object of class TS

# Regenerate the data in TS under the null distribution 
setGeneric("H0_gen", function(object) {
  standardGeneric("H0_gen")
})


setMethod("H0_gen",signature(object = "TS"), function(object) {
  stop("The H0_gen function generate data under the null hypothesis that
       depends on the type of analysis. Use or implement and appropriate class")
})

#' Generate data under an alternative distribution
#' @description generic function to generate data under an alternative distribution. The generation starts from the null and modifies generation parameters for the leaves in \code{alternativeLeaves}. The argument \code{alternativeParameter} contains parameters used to change the distribution of the \code{alternativeLeaves}. The meaning of this parameter varies according to the specific TBSS. Instances of this function are used to compute power in \code{evalPower}. 
#' @param object an object of class TS
#' @param alternativeParameter a number or a vector of length \code{length(alternativeLeaves)} containing parameters to be used to change the distribution of the \code{alternativeLeaves}
#' @param alternativeLeaves a vector containing the IDs of the leaves to be modified from the null distribution.

setGeneric("H1_gen", function(object, alternativeParameter, alternativeLeaves) {
  standardGeneric("H1_gen")
})

setMethod("H1_gen",signature(object = "TS"), function(object, alternativeParameter, alternativeLeaves){
  stop("The H1_gen function generate data under the null hypothesis that
       depends on the type of analysis. Use or implement and appropriate class")
})


# ---deprecated! to remove soon
setGeneric("bootstrap", function(object,parallel = FALSE, ncpus = NULL) {
  standardGeneric("bootstrap")
})

#setMethod("bootstrap", signature(object = "TS"), function(object) {
#  stop("The bootstrap algorithm depends on the type of analysis. Use or implement and appropriate class")
#})

#' Monte Carlo algorithm 
#' @description Function that approximates the null distribution of the node-specific discrepancy measures. Results are stored in the slot \code{LRT_H0} 
#' @param object an object of class TS
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.

  monteCarlo = function(object, parallel = FALSE, ncpus = NULL)
  {
   # a single MC sample
   do_MC = function(...)
   {
    boot_object = object
    ## re-generate data under the null
    boot_object = H0_gen(object)
 
    retval = computeLRT(boot_object)@LRT
    return(retval)
}
if(parallel)
{
 suppressMessages(require(parallel))
 suppressMessages(require(TBSS))
 cl <- makeCluster(ncpus)
 clusterExport(cl, "do_MC", envir = environment())
 clusterExport(cl, "object", envir = environment())
 boot_LRT = parSapply(cl, 1:object@B, do_MC)
 stopCluster(cl)
}
else
{
	## non-parallel version
 boot_LRT =  sapply(1:object@B, do_MC)
}

 object@LRT_H0 = boot_LRT
 return(object)
}
  


setMethod(f = "show",
          signature = "TS",
          definition = function(object){

          ranking = order(object@LRT,decreasing = TRUE)
          boot_maxs = apply(object@LRT_H0,2,max, na.rm = TRUE)
          pvalues = sapply(object@LRT[ranking], \(x) (sum(boot_maxs >= x) +1)/(object@B +1))
          #### most likely cut  
          sprintf("The most likely cut is %s, with log-likelihood ratio = %.3f, and p-value = %.3f",       names(object@LRT)[ranking[1]], max(object@LRT),  pvalues[1]) |> print()
          })



#' @export
setMethod(f = "summary",
          signature = "TS",
          definition = function(object){

          ranking = order(object@LRT,decreasing = TRUE)
          boot_maxs = apply(object@LRT_H0,2,max)
          pvalues = sapply(object@LRT[ranking], \(x) sum(boot_maxs >= x)/(object@B +1))

	  
	  retval = data.frame(node = object@nodesToTest[ranking],
		                     LRT = object@LRT[ranking],
		                     pvalue =  pvalues)
	 
          rownames(retval) = NULL
	  return(retval)
	  })

