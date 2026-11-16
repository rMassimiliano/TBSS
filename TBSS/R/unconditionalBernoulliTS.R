#############################################
#############################################
### Unconditional Bernoulli Tree Scan
#############################################
#############################################

#' Sufficient Statistics
#' Compute sufficient statistics for unconditionalBernoulliTS
#' @param object an object of class TS
#' @rdname unconditionalBernoulliTS

setMethod("computeSS", signature(object = "unconditionalBernoulliTS"), function(object)
{
  object@nodeSS = data.frame(node = names(object@mapNodesLeaves), case = 0L, control = 0L)

  for(r in 1:length(object@mapNodesLeaves))
  {
   ind = which(object@leaves %in% object@mapNodesLeaves[[r]])
   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   object@nodeSS$case[r]     =  sum(object@data$case[ind])
   object@nodeSS$control[r]  = sum(object@data$control[ind])
   
  }
  return(object)
}
)

#' log-likelihood ratio statistics 
#' Compute The log-likelihood ratio statistics for unconditionalBernoulliTS
#' @param object an object of class TS
#' @rdname unconditionalBernoulliTS
setMethod("computeLRT", signature(object = "unconditionalBernoulliTS"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## compute LLR for all the nodes
  LRT = numeric(length(object@nodesToTest))
  j_inds = match(object@nodesToTest, object@nodeSS$node) ## row of nodeSS for each tested node, computed outide the loop
  for(r in 1:length(object@nodesToTest))
  {
   j = j_inds[r]
   LRT[r] =  uncon_bernoulli_llrt(object@nodeSS$case[j],object@nodeSS$control[j],object@p)
  }
   names(LRT) = object@nodesToTest
   object@LRT = LRT
   return(object)
 })


#' generate a new tree under the global null H0
#' @param object an object of class TS
#' @rdname unconditionalBernoulliTS
setMethod("H0_gen", signature(object = "unconditionalBernoulliTS"),
function(object)
{
    ## generate the leaves data under H0
    ## we generate data from the leaves even if these are not tested	  
     leaves_list =  object@leaves[which(object@leaves %in% object@nodesToTest)]
     leaves_list = object@leaves   

     ## update leaves
     indx = which(object@nodeSS$node %in% leaves_list)
     totals = object@nodeSS$case[indx] + object@nodeSS$control[indx] 
     cases  = stats::rbinom(length(indx), size = totals, prob = object@p)
     object@nodeSS$case[indx] = cases
     object@nodeSS$control[indx] = totals - cases

     # update other nodes
     ind_set = which(!(object@nodeSS$node%in% leaves_list))
     ind_set = ind_set[which(object@nodeSS$node[ind_set] %in% object@nodesToTest)]
     for(r in ind_set)
     {
      leaf = object@mapNodesLeaves[[object@nodeSS$node[r]]]
      object@nodeSS$case[r] =  sum(object@nodeSS$case[object@nodeSS$node %in% leaf])
      object@nodeSS$control[r] =  sum(object@nodeSS$control[object@nodeSS$node %in% leaf])
    }
    return(object)
  }
)



#'  Summary method
#' @description provide a  summary of the results of and unconditional Bernoulli TBSS analysis. 
#' @param object an object of class TS
#' @rdname unconditionalBernoulliTS
#' @importFrom rlang .data
setMethod(f = "summary",
          signature = "unconditionalBernoulliTS",
          definition = function(object){

          ranking = order(object@LRT,decreasing = TRUE)
          boot_maxs = apply(object@LRT_H0,2,max)
          pvalues = sapply(object@LRT[ranking], \(x) (sum(boot_maxs >= x)+1)/(object@B +1))

	  retval = data.frame(node = names(object@LRT)[ranking],
		                     LRT = object@LRT[ranking],
		                     pvalue =  pvalues)

          suppressMessages({
		  retval = retval |> dplyr::full_join(object@nodeSS);
                  retval = retval |>
	           dplyr::mutate(expected = (.data$control + .data$case)*object@p) |>
		   dplyr::mutate(relativeRisk =  .data$case/.data$expected, excessCases = .data$case - .data$expected) |>
		   dplyr::select(.data$node,.data$case, .data$control, .data$expected, .data$relativeRisk, .data$excessCases, .data$LRT,.data$pvalue);
		  retval = retval |> 
			  dplyr::arrange(.data$pvalue);
		  rownames(retval) = NULL;
	          })

	  return(retval)
	  })

