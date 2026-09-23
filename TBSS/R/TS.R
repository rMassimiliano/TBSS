#' Populate the leaves slot
#' @description Populate the slot leaves of TS with a vector of the unique values of the first column of data (leaves ID). 
#' @param object an object of class TS
#' @rdname TS
listLeaves = function(object)
{
 object@leaves =  unique(object@data[,1])
 return(object)
}

#' Populate the slot nodes
#' @description Populate the slot node of TS with a vector containing the unique value of the first column of the slot \code{tree}. Empty \code{character} are not allowed
#' @param object an object of class TS
#' @rdname TS
listNodes = function(object)
{
  object@nodes =  setdiff(unique(unlist(object@tree)),"")
  names(object@nodes) = NULL
  object
}


#' Map nodes to the corresponding leaves 
#' @description Create a named list. Each element of the list is named after a node and contains the ID of leaves linked to that node. The function returns an object of class TS with a populated slot \code{mapNodeLeaves}.
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
#' @description generic function to compute node-specific sufficient statistics. The function populates the slot \code{nodeSS}. 
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
#' @description generic function to compute node-specific discrepancy measures. The function populates the slot \code{nodeSS}. 
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
#' @description generic function to generate data under the null distribution. Instances of this function are used in \code{monteCarlo} 
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
#' @description generic function to generate data under an alternative distribution. The generation starts from the null and modifies generation parameters for the leaves in \code{alternativeLeaves}. The argument \code{alternativeParameter} contains parameters used to change the distribution of the \code{alternativeLeaves}. The meaning of this parameter varies according to the specific TBSS. Instances of this function are used to compute power in \code{evalPower}. 
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
 parallel::clusterExport(cl, "do_MC", envir = environment())
 parallel::clusterExport(cl, "object", envir = environment())
 boot_LRT = parallel::parSapply(cl, 1:object@B, do_MC)
 parallel::stopCluster(cl)
}
else
{
	## non-parallel version
 boot_LRT =  sapply(1:object@B, do_MC)
}

 object@LRT_H0 = boot_LRT
 return(object)
}
  


#'  Show method
#' @description  Print the most likely cut in a TBSS analysis with log-likelihood ratio and p-value
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
#' @description provide a (minimal) summary of the TBSS analysis results.  Specific instances of TS class might implement more detailed summaries
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

