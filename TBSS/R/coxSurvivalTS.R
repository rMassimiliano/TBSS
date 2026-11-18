#' @importFrom rlang .data
#' @importFrom survival Surv coxph
setClass(Class = "coxSurvivalTS",  contains = 'permutationPoissonTS')

#' Sufficient Statistics
#' Compute sufficient statistics for the permutationPoissonTS 
#' @param object an object of class TS
#' @rdname permutationPoissonTS
setMethod("computeSS", signature(object = "coxSurvivalTS"), function(object)
{
  object@nodeSS = data.frame(node = names(object@mapNodesLeaves), observed = 0L, observed_pt = 1L, expected = 0L, expected_pt = 1L, HR = 0L, LRT = 0L)

  for(r in 1:length(object@mapNodesLeaves))
  {
   d = object@idInNode[[names(object@mapNodesLeaves)[r]]]
   p = object@personTime[[names(object@mapNodesLeaves)[r]]]
   q = object@personTimeNotInNode[[names(object@mapNodesLeaves)[r]]]

   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   object@nodeSS$observed[r] = sum(object@idInNode[[r]]$weight*(object@idInNode[[r]]$exposure == 1))
   object@nodeSS$expected[r] = sum( object@idInNode[[r]]$weight* (object@idInNode[[r]]$exposure == 0))

   object@nodeSS$expected_pt[r] = sum((d$weight*d$followup)[d$exposure == 0]) + sum((p$weight*p$followup)[p$exposure == 0]) + sum((q$weight*q$followup)[q$exposure == 0])

   object@nodeSS$observed_pt[r] = sum((d$weight*d$followup)[d$exposure == 1]) + sum((p$weight*p$followup)[p$exposure == 1]) + sum((q$weight*q$followup)[q$exposure == 1])

    ## estimate cox model

    tmpData  = data.frame(time = c(p$followup, d$followup, q$followup),
                          exposure = c(p$exposure, d$exposure, q$exposure),
                          weight = c(p$weight, d$weight, q$weight),
                          status = rep(c(0, 1, 0), c(NROW(p), NROW(d), NROW(q))))

   
    ## estimate cox model
  
  	if(NROW(tmpData)>2)
  	{
  		cmod = survival::coxph(survival::Surv(time,status)~exposure, weight = tmpData$weight, data = tmpData, ties = "breslow")
		object@nodeSS$HR[r] = exp(cmod$coef[[1]])
		object@nodeSS$LRT[r] = -2 * (cmod$loglik[1] - cmod$loglik[2])
  
  	}
	else
	{
		object@nodeSS$HR[r] = NA
		object@nodeSS$LRT[r] = NA
	}

  }

  
  return(object)
}
)


#' log-likelihood ratio statistics 
#' Compute The log-likelihood ratio statistics for permutationPoissonTS
#' @param object an object of class TS
#' @rdname permutationPoissonTS
setMethod("computeLRT", signature(object = "coxSurvivalTS"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## this is pre-compute in the sufficient statistic node. 
  L = length(object@nodesToTest)
  LRT = numeric(L)
  coefs = numeric(L)
  i=1
  for(node in object@nodesToTest)
  {
   j = match(node,object@nodeSS$node)
   LRT[i] = object@nodeSS$LRT[j]
   coefs[i] = log(object@nodeSS$HR[j])
   i=i+1
  }

   LRT = applyDirection(object, LRT, coefs)
   object@LRT = LRT
   names(object@LRT) = object@nodesToTest
   return(object)
 })




#'  Summary method
#' @description provide a  summary of the TBSS analysis results. Include Hazard Ratio (HR), time and relative measures 
#' @param object an object of class TS

#' @return \code{summary()} returns a \code{data.frame}. One row per node of the tree: the tested nodes first, ordered by decreasing \code{LRT}, then the nodes that were not tested, with \code{LRT} and \code{pvalue} equal to \code{NA}.  \code{observed_rate} and \code{expected_rate} are the weighted event rates (events per unit of person time) of the exposed and of the unexposed; \code{relativeRisk} is their ratio and \code{excessCases} their difference. \code{pvalue} is the Monte Carlo p-value adjusted for the multiple testing over the tree, (1 + number of replicates whose maximum LRT is at least the LRT of the node)/(B + 1). \code{HR}, is the hazard ratio of the exposure estimated by the weighted Cox model of the node.
#' @rdname permutationPoissonTS
setMethod(f = "summary",
          signature = "coxSurvivalTS",
          definition = function(object){

          ranking   = order(object@LRT,decreasing = TRUE)
          boot_maxs = apply(object@LRT_H0,2,max, na.rm = TRUE)
          pvalues   = sapply(object@LRT[ranking], \(x) (sum(boot_maxs >= x) +1)/(object@B +1))


	  retval = data.frame(node = names(object@LRT)[ranking],
		                     LRT = object@LRT[ranking],
		                     pvalue =  pvalues)
          rownames(retval) = NULL
	   suppressMessages({
	retval =   retval |>
	           dplyr::full_join(object@nodeSS) |>
	            dplyr::mutate(HR= .data$HR, relativeRisk = (.data$observed/.data$observed_pt)/(.data$expected/.data$expected_pt), excessCases = .data$observed/.data$observed_pt - .data$expected/.data$expected_pt) |>
		    dplyr::mutate(observed_rate = .data$observed/.data$observed_pt, expected_rate = .data$expected/.data$expected_pt) |>
	            dplyr::select(.data$node, .data$HR, .data$observed_rate, .data$expected_rate, .data$relativeRisk, .data$excessCases,.data$LRT,.data$pvalue) 
	  })
	  return(retval)
	  })
