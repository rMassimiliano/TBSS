#############################################
#############################################
### Conditional Possion Tree Scan
#############################################
#############################################
#' @importFrom rlang .data
setMethod("computeLRT", signature(object = "conditionalPoissonTS"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## compute LLR for all the nodes
  TOTAL_EXPECTED = sum(object@data$expected)
  TOTAL_OBSERVED = sum(object@data$observed)

  LRT = numeric(length(object@nodesToTest))
  for(r in 1:length(object@nodesToTest))
  {
   j = which(object@nodeSS$node == object@nodesToTest[r])
   LRT[r] = uncon_poisson_llrt(object@nodeSS$observed[j],object@nodeSS$expected[j] * TOTAL_OBSERVED/TOTAL_EXPECTED)
  }
   object@LRT = LRT
   names(object@LRT) = object@nodesToTest
   return(object)
 })



##setMethod("computeSS", signature(object = "conditionalPoissonTS"), function(object)
##{
##  object@nodeSS = data.frame(node = names(object@mapNodesLeaves), observed = 0L, expected = 0L)
##
##  TOTAL_EXPECTED = sum(object@data$expected)
##  TOTAL_OBSERVED = sum(object@data$observed)
##  
##  for(r in 1:length(object@mapNodesLeaves))
##  {
##   ind = which(object@leaves %in% object@mapNodesLeaves[[r]])
##   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
##   object@nodeSS$observed[r] = sum(object@data$observed[ind])
##   object@nodeSS$expected[r] = sum(object@data$expected[ind]) ##*TOTAL_OBSERVED/TOTAL_EXPECTED
##  }
##  return(object)
##}
##)







#' generate a new tree under the global null H0
#' @param object an object of class TS
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


#' Summary method
#' @description provide a (minimal) summary of the TBSS analysis results.  Specific instances of TS class might implement more detailed summaries
#' @param object an object of class TS
#'@rdname TS
#'@export
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
