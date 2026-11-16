## This file pre-define the classes used in the package to avoid issues due to the fact that files are loaded in alphabetical order
#General TS class for tree scan
#'General TS class for TBSS
#' General TS class.
#' @description This class is used as a template for all the TBSS implementation. Several generic function (e.g. utilities to parse the tree) are implemented in this generic class
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
#'@slot direction a \code{character} define the direction of the test. For one-side test LRT in the side that is not considered is set to 0 both for the observed data and the ones replicated under the null.
#' @docType class
#' @name TS-class
#' @aliases TS-class TS
setClass(Class = "TS", slots = c(data = 'data.frame',
                                 tree ='data.frame',
                                 leaves = 'character',
                                 nodes = 'character',
                                 mapNodesLeaves = 'list',
				 nodesToTest ='character',
                                 nodeSS = 'data.frame',
                                 LRT ='numeric',
                                 B = 'numeric',
                                 LRT_H0 = 'matrix',
				 direction = 'character'))
#' Utility functions of the TS class
#' @description \code{listLeaves()} Populate the slot leaves of TS with a vector of the unique values of the first column of data (leaves ID). 
#' @param object an object of class TS
#' @rdname TS
listLeaves = function(object)
{
 if(anyNA(object@data[,1])) stop("The first column of the data (leaf) contains missing values.")
 object@leaves =  unique(object@data[,1])
 return(object)
}


#' Populate the slot nodes
#' @description \code{listNodes()} Populate the slot node of TS with a vector containing the unique value of the first column of the slot \code{tree}. Empty \code{character} are not allowed
#' @param object an object of class TS
#' @rdname TS
listNodes = function(object)
{
  object@nodes =  setdiff(unique(unlist(object@tree)),"")
  names(object@nodes) = NULL
  object
}


#' Map nodes to the corresponding leaves 
#' @description \code{mapNodesLeaves()} Create a named list. Each element of the list is named after a node and contains the ID of leaves linked to that node. The function returns an object of class TS with a populated slot \code{mapNodeLeaves}.
#' @param object an object of class TS
#' @param ...  additional arguments (currently for retro-compatibility, not used)
#' @rdname TS
mapNodesLeaves =  function(object, ...)
{
  if(length(object@leaves) == 0) object = listLeaves(object)
  if(length(object@nodes) == 0) object = listNodes(object)
	
  nodeList =  list()
  # add leaves as themselves
  for(i in seq_along(object@leaves))
  {
    nodeList[[object@leaves[i]]] = object@leaves[i]
  }
  ## initialized explored nodes (bfs) and not visited nodes
  bfs = object@leaves
  notvisited = setdiff(object@nodes, object@leaves) ## at starting point we only visited the leaf nodes

  # start the tree exploration
  while(TRUE)
  {
   w          = object@tree[,1] %in% bfs
   nextlevel  = unique(object@tree[w,2])
   nextlevel  = intersect(nextlevel,notvisited)
   notvisited = setdiff(notvisited,nextlevel)

   if(all(nextlevel=="")) break
   # go through all edges in w and add corresponding leaves
		for(i in which(w)) {
			child              = object@tree[i,1]
			parent             = object@tree[i,2]
			temp               = c(nodeList[[parent]],nodeList[[child]])
			nodeList[[parent]] = unique(temp)
		}
		# next bfs sweep
		bfs = nextlevel
   }
 object@mapNodesLeaves = nodeList
 node_removed  = length(unique(unlist(object@tree))) - length(nodeList)
  if(node_removed>0) warning(sprintf("We removed %i nodes from the analysis because no data have been provided for the corresponding leaves", node_removed))
  return(object) 
}



#' Compute node-specific sufficient statistics
#' @description \code{computeSS()} generic function to compute node-specific sufficient statistics. The function populates the slot \code{nodeSS}. 
#' @param object an object of class TS
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.
#' @rdname TS
setGeneric("computeSS", function(object, parallel = FALSE, ncpus = NULL) {
  standardGeneric("computeSS")
})

setMethod("computeSS", signature(object = "TS"), function(object,parallel = FALSE, ncpus = NULL) {
  stop("To compute node specific sufficient statistics we need a model. Use or implement an appropriate function.")
})




#' Compute node-specific discrepancy measures
#' @description \code{computeLRT()} generic function to compute node-specific discrepancy measures. The function populates the slot \code{nodeSS}. 
#' @param object an object of class TS
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.
#' @rdname TS
setGeneric("computeLRT", function(object, parallel = FALSE, ncpus = NULL) {
  standardGeneric("computeLRT")
})

# Compute the LRT. The expression depends on the model
setMethod("computeLRT", signature(object = "TS"), function(object,parallel = FALSE, ncpus = NULL) {
  stop("The LRT formula depends on the type of analysis. Use or implement an appropriate class.")
})

#' Generate data under the null distribution
#' @description \code{H0_gen()} generic function to generate data under the null distribution. Instances of this function are used in \code{monteCarlo} 
#' @param object an object of class TS
# Regenerate the data in TS under the null distribution 
#' @rdname TS
setGeneric("H0_gen", function(object) {
  standardGeneric("H0_gen")
})


setMethod("H0_gen",signature(object = "TS"), function(object) {
  stop("The H0_gen function generate data under the null hypothesis that
       depends on the type of analysis. Use or implement and appropriate class")
})

#' Generate data under an alternative distribution
#' @description \code{H1_gen()} generic function to generate data under an alternative distribution. The generation starts from the null and modifies generation parameters for the leaves in \code{alternativeLeaves}. The argument \code{alternativeParameter} contains parameters used to change the distribution of the \code{alternativeLeaves}. The meaning of this parameter varies according to the specific TBSS. Instances of this function are used to compute power in \code{evalPower}. 
#' @param object an object of class TS
#' @param alternativeParameter a number or a vector of length \code{length(alternativeLeaves)} containing parameters to be used to change the distribution of the \code{alternativeLeaves}
#' @param alternativeLeaves a vector containing the IDs of the leaves to be modified from the null distribution.
#' @rdname TS
setGeneric("H1_gen", function(object, alternativeParameter, alternativeLeaves) {
  standardGeneric("H1_gen")
})

setMethod("H1_gen",signature(object = "TS"), function(object, alternativeParameter, alternativeLeaves){
  stop("The H1_gen function generate data under the null hypothesis that
       depends on the type of analysis. Use or implement and appropriate class")
})


#' Monte Carlo algorithm 
#' @description Function that approximates the null distribution of the node-specific discrepancy measures. Results are stored in the slot \code{LRT_H0} 
#' @param object an object of class TS
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.
#' @rdname TS
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
 cl = parallel::makeCluster(ncpus)
 on.exit(parallel::stopCluster(cl), add = TRUE)
 parallel::clusterSetRNGStream(cl, sample.int(.Machine$integer.max, 1))
 parallel::clusterExport(cl, "do_MC", envir = environment())
 parallel::clusterExport(cl, "object", envir = environment())
 boot_LRT = parallel::parSapply(cl, 1:object@B, do_MC)
# parallel::stopCluster(cl)
}
else
{
	## non-parallel version
 boot_LRT =  sapply(1:object@B, do_MC)
}
 if(!is.matrix(boot_LRT)) boot_LRT = matrix(boot_LRT, nrow = 1) ## exception for testing one single node 
 object@LRT_H0 = boot_LRT
 return(object)
}
  


#'  Show method
#' @description  \code{show()} print the most likely cut in a TBSS analysis with log-likelihood ratio and p-value
#' @param object an object of class TS
#' @rdname TS
setMethod(f = "show",
          signature = "TS",
          definition = function(object){

          ranking = order(object@LRT,decreasing = TRUE)
          boot_maxs = apply(object@LRT_H0,2,max, na.rm = TRUE)
          pvalues = sapply(object@LRT[ranking], \(x) (sum(boot_maxs >= x) +1)/(object@B +1))
          #### most likely cut  
          sprintf("The most likely cut is %s, with log-likelihood ratio = %.3f, and p-value = %.3f",       names(object@LRT)[ranking[1]], max(object@LRT),  pvalues[1]) |> print()
          })


#'  Summary method
#' @description \code{summary()} provide a (minimal) summary of the TBSS analysis results.  Specific instances of TS class might implement more detailed summaries
#' @param object an object of class TS
#' @rdname TS
#'@export
setMethod(f = "summary",
          signature = "TS",
          definition = function(object){

          ranking = order(object@LRT,decreasing = TRUE)
          boot_maxs = apply(object@LRT_H0,2,max)
          pvalues = sapply(object@LRT[ranking], \(x) (sum(boot_maxs >= x) +1)/(object@B +1))


	  
	  retval = data.frame(node = object@nodesToTest[ranking],
		                     LRT = object@LRT[ranking],
		                     pvalue =  pvalues)
	 
          rownames(retval) = NULL
	  return(retval)
	  })




#' One-sided restriction of the discrepancy measures
#' @description \code{applyDirection()} for one sided tests sets to 0 the LRT in the direction that is not considered in the test. For \code{direction == 'positive'} the LRT with \code{excess<=0} are set to 0 while for \code{direction == 'negative'} the LRT with \code{excess>=0} are set to 0. The function is intended to be used within \code{computeLRT()} so the same restriction applies to the observed data and to every Monte Carlo replicate.
#' @param object an object of class TS
#' @param LRT a vector of node-specific LRT 
#' @param excess a vector, of the same length as \code{LRT}, with the node-specific excess (e.g., observed minus expected)
#' @rdname TS
applyDirection = function(object, LRT, excess)
{
  if(length(object@direction) == 0) return(LRT)
  if(object@direction == 'positive') LRT[excess <= 0] = 0
  else if(object@direction == 'negative') LRT[excess >=0] = 0
  return(LRT)
}
