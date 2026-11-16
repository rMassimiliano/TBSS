#############################################
#############################################
### Conditional Possion Tree Scan
#############################################
#############################################
#' @importFrom rlang .data
setClass(Class = "conditionalPoissonTS",  contains = 'unconditionalPoissonTS')

#' log-likelihood ratio statistics 
#' @description \code{computeLRT()} Compute The log-likelihood ratio statistics for conditionalPoissonTS 
#' @param object an object of class TS
#' @rdname conditionalPoissonTS 
setMethod("computeLRT", signature(object = "conditionalPoissonTS"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## compute LLR for all the nodes
  leaves_ind = object@nodeSS$node %in% object@leaves
  TOTAL_EXPECTED = sum(object@nodeSS$expected[leaves_ind])
  TOTAL_OBSERVED = sum(object@nodeSS$observed[leaves_ind])

  LRT = numeric(length(object@nodesToTest))
  excess = numeric(length(object@nodesToTest))
  j_inds = match(object@nodesToTest, object@nodeSS$node) ## row of nodeSS for each tested node, computed outide the loop
  for(r in 1:length(object@nodesToTest))
  {
   j = j_inds[r]
   LRT[r] = uncon_poisson_llrt(object@nodeSS$observed[j],object@nodeSS$expected[j] * TOTAL_OBSERVED/TOTAL_EXPECTED)
   excess[r] = object@nodeSS$observed[j] - object@nodeSS$expected[j] * TOTAL_OBSERVED/TOTAL_EXPECTED
  }
   LRT = applyDirection(object, LRT, excess)
   object@LRT = LRT
   names(object@LRT) = object@nodesToTest
   return(object)
 })


#' Data under the null model
#' @description \code{H0_gen()} generate data under the null
#' @param object an object of class TS
#' @rdname conditionalPoissonTS 
setMethod("H0_gen", signature(object = "conditionalPoissonTS"),
function(object)
{

  TOTAL_EXPECTED = sum(object@data$expected)
  TOTAL_OBSERVED = sum(object@data$observed)
  ## generate the leaves data under H0
  leaves_list = object@leaves   

  indx = which(object@nodeSS$node %in% leaves_list)
  probs =  object@nodeSS$expected[indx]/sum(object@nodeSS$expected[indx])

  object@nodeSS$observed[indx] =    stats::rmultinom(1,size = TOTAL_OBSERVED, prob = probs)
  
     ## update only the SS for the nodes that are used in the LRT
  ind_set = which(!(object@nodeSS$node%in% leaves_list))
  ind_set = ind_set[which(object@nodeSS$node[ind_set] %in% object@nodesToTest)]
   for(r in ind_set)
   {
    leaf = object@mapNodesLeaves[[object@nodeSS$node[r]]]
    object@nodeSS$observed[r] =  sum(object@nodeSS$observed[object@nodeSS$node %in% leaf])
      object@nodeSS$expected[r] =  sum(object@nodeSS$expected[object@nodeSS$node %in% leaf])
    }
    return(object)
}
)
#' Data under an alternative model
#' @description \code{H1_gen()} generate data under a specific alternative model, conditional on the total number of observed events.  The expected counts of the \code{alternativeLeaves} are multiplied by \code{alternativeParameter}. Events are generated from  a multinomial distribution.
#' @param object an object of class TS
#' @param alternativeParameter relative risk(s) used for the leaves under the alternative
#' @param alternativeLeaves leaves under the alternative model
#' @rdname conditionalPoissonTS 
setMethod("H1_gen", signature(object = "conditionalPoissonTS"),
function(object, alternativeParameter, alternativeLeaves)
{

  TOTAL_OBSERVED = sum(object@data$observed)
  ## generate the leaves data under H1
  leaves_list = object@leaves   

  indx = which(object@nodeSS$node %in% leaves_list)
  ## relative risk of each leaf: 1 under the null, alternativeParameter for the alternative leaves
  rr = rep(1, length(indx))
  alt = match(alternativeLeaves, object@nodeSS$node[indx])
  rr[alt] = rep_len(alternativeParameter, length(alt))
  probs =  object@nodeSS$expected[indx]*rr/sum(object@nodeSS$expected[indx]*rr)

  object@nodeSS$observed[indx] =    stats::rmultinom(1,size = TOTAL_OBSERVED, prob = probs)
  
     ## update only the SS for the nodes that are used in the LRT
  ind_set = which(!(object@nodeSS$node%in% leaves_list))
  ind_set = ind_set[which(object@nodeSS$node[ind_set] %in% object@nodesToTest)]
   for(r in ind_set)
   {
    leaf = object@mapNodesLeaves[[object@nodeSS$node[r]]]
    object@nodeSS$observed[r] =  sum(object@nodeSS$observed[object@nodeSS$node %in% leaf])
      object@nodeSS$expected[r] =  sum(object@nodeSS$expected[object@nodeSS$node %in% leaf])
    }
    return(object)
}
)


#' Summary method
#' @description \code{summary()} provide a summary of the conditional poisson TBSS analysis results.  
#' @param object an object of class TS
#' @rdname conditionalPoissonTS 
#' @export 
setMethod(f = "summary",
          signature = "conditionalPoissonTS",
          definition = function(object){
          TOTAL_EXPECTED = sum(object@data$expected)
          TOTAL_OBSERVED = sum(object@data$observed)

          ranking = order(object@LRT,decreasing = TRUE)
          boot_maxs = apply(object@LRT_H0,2,max)
          pvalues   = sapply(object@LRT[ranking], \(x) (sum(boot_maxs >= x) +1)/(object@B +1))


	  retval = data.frame(node = names(object@LRT)[ranking],
		                     LRT = object@LRT[ranking],
		                     pvalue =  pvalues)
          rownames(retval) = NULL
	  suppressMessages({
	retval =   retval |>
	           dplyr::full_join(object@nodeSS) |>
	           dplyr::mutate(relativeRisk = (.data$observed/.data$expected)/((TOTAL_OBSERVED - .data$observed)/(TOTAL_EXPECTED - .data$expected)), excessCases = .data$observed - .data$expected *(TOTAL_OBSERVED -.data$observed)/(TOTAL_EXPECTED - .data$expected)) |>
		   dplyr::mutate(expected = .data$expected * TOTAL_OBSERVED/TOTAL_EXPECTED) |>
	           dplyr::select(.data$node,.data$observed,.data$expected,.data$relativeRisk, .data$excessCases,.data$LRT,.data$pvalue) 
	  })
	  return(retval)
	  })
