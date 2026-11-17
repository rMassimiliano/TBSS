#############################################
#############################################
### Permutation Possion Tree Scan
#############################################
#############################################
#' @importFrom rlang .data
setClass(Class = "permutationPoissonTS",  contains = 'TS', slots = c(idInNode = 'list',personTime = 'list',personTimeNotInNode = 'list', cohort = 'data.frame', fup = 'logical')) 

#' Sufficient Statistics
#' Compute sufficient statistics for the permutationPoissonTS 
#' @param object an object of class TS
#' @rdname permutationPoissonTS
setMethod("computeSS", signature(object = "permutationPoissonTS"), function(object)
{
  object@nodeSS = data.frame(node = names(object@mapNodesLeaves), observed = 0L, observed_pt = 1L, expected = 0L, expected_pt = 1L)

  for(r in 1:length(object@mapNodesLeaves))
  {
   d = object@idInNode[[names(object@mapNodesLeaves)[r]]]
   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   object@nodeSS$observed[r] = sum(d$weight*(d$exposure == 1))
   object@nodeSS$expected[r] = sum( d$weight*(d$exposure == 0))
   if(object@fup)
   {
    ## combine person time for censored and not in node
     p = object@personTime[[names(object@mapNodesLeaves)[r]]] 
     q = object@personTimeNotInNode[[names(object@mapNodesLeaves)[r]]]
     object@nodeSS$expected_pt[r] = sum((d$weight*d$followup)[d$exposure == 0]) + sum((p$weight*p$followup)[p$exposure == 0]) + sum((q$weight*q$followup)[q$exposure == 0])
    object@nodeSS$observed_pt[r] = sum((d$weight*d$followup)[d$exposure == 1]) + sum((p$weight*p$followup)[p$exposure == 1]) + sum((q$weight*q$followup)[q$exposure == 1])
   }

  }

  
  return(object)
}
)


#' log-likelihood ratio statistics 
#' Compute The log-likelihood ratio statistics for permutationPoissonTS
#' @param object an object of class TS
#' @rdname permutationPoissonTS
setMethod("computeLRT", signature(object = "permutationPoissonTS"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## compute LLR for all the nodes
  LRT = numeric(length(object@nodesToTest))
  excess = numeric(length(object@nodesToTest))
  j_inds = match(object@nodesToTest, object@nodeSS$node) ## row of nodeSS for each tested node, computed outide the loop
  for(r in 1:length(object@nodesToTest))
  {
   j = j_inds[r]
   LRT[r] = with(object@nodeSS[j,], poisson_llrt(observed, observed_pt, expected, expected_pt))
   excess[r] = with(object@nodeSS[j,], observed/observed_pt - expected/expected_pt)
  }
   LRT = applyDirection(object, LRT, excess)
   object@LRT = LRT
   names(object@LRT) = object@nodesToTest
   return(object)
 })

#' generate data under the null via permutation of the lables
#' @param object an object of class TS
#' @rdname permutationPoissonTS
setMethod("H0_gen", signature(object = "permutationPoissonTS"),
function (object)
{
    ## generate the leaves data under H0

	## it is ordered by strata
  #permuted_cohort = object@cohort |> dplyr::select(.data$patientID,.data$strata) |>  dplyr::group_by(.data$strata) |> dplyr::slice_sample(prop = 1)

 ## cohort is ordered by strata so this change allocation but not structure since the PatientId are permuted withing the strata
 #permuted_cohort$exposure = object@cohort$exposure
 #permuted_cohort$weight = object@cohort$weight
 
  perm = seq_len(NROW(object@cohort))
  perm[order(object@cohort$strata)] = order(object@cohort$strata, stats::runif(NROW(object@cohort))) ## use a random order within strata
  permuted_cohort = data.frame(patientID = object@cohort$patientID[perm], exposure = object@cohort$exposure, weight = object@cohort$weight)

     ## update node that have to be tested
     for(j in 1:length(object@nodesToTest))
     {
	r = object@nodesToTest[j] 
     ## change exposure and weights for the SS order is given by nodes in idNode
 	m = match(object@idInNode[[r]]$patientID, permuted_cohort$patientID)
        object@idInNode[[r]]$exposure = permuted_cohort$exposure[m]
        object@idInNode[[r]]$weight   = permuted_cohort$weight[m]
      if(object@fup)
      {
        m = match(object@personTime[[r]]$patientID, permuted_cohort$patientID)
        object@personTime[[r]]$exposure = permuted_cohort$exposure[m]
        object@personTime[[r]]$weight   = permuted_cohort$weight[m]

        m = match(object@personTimeNotInNode[[r]]$patientID, permuted_cohort$patientID)
        object@personTimeNotInNode[[r]]$exposure = permuted_cohort$exposure[m]
        object@personTimeNotInNode[[r]]$weight   = permuted_cohort$weight[m]
      }
     }
     ## compute SS  --- only for nodes to be tested
     object@mapNodesLeaves = object@mapNodesLeaves[object@nodesToTest]
     object@idInNode = object@idInNode[object@nodesToTest]
     if(object@fup)
     {
      object@personTime = object@personTime[object@nodesToTest]
      object@personTimeNotInNode = object@personTimeNotInNode[object@nodesToTest]
     }
     object = computeSS(object)
    return(object)
   }
)
	  

#'  Summary method
#' @description provide a  summary of the TBSS analysis results. Include time and relative measure if follow-up time is used. 
#' @param object an object of class TS
#' @rdname permutationPoissonTS

setMethod(f = "summary",
          signature = "permutationPoissonTS",
          definition = function(object){

          ranking   = order(object@LRT,decreasing = TRUE)
          boot_maxs = apply(object@LRT_H0,2,max, na.rm = TRUE)
          pvalues   = sapply(object@LRT[ranking], \(x) (sum(boot_maxs >= x) +1)/(object@B +1))


	  retval = data.frame(node = names(object@LRT)[ranking],
		                     LRT = object@LRT[ranking],
		                     pvalue =  pvalues)
          rownames(retval) = NULL
	  if(object@fup)
	  {
	   suppressMessages({
	retval =   retval |>
	           dplyr::full_join(object@nodeSS) |>
	            dplyr::mutate(relativeRisk = (.data$observed/.data$observed_pt)/(.data$expected/.data$expected_pt), excessCases = .data$observed/.data$observed_pt - .data$expected/.data$expected_pt) |>
		    dplyr::mutate(observed_rate = .data$observed/.data$observed_pt, expected_rate = .data$expected/.data$expected_pt) |>
	            dplyr::select(.data$node,.data$observed_rate, .data$expected_rate, .data$relativeRisk, .data$excessCases,.data$LRT,.data$pvalue) 
	  })
	  }
	  else
	  {
	   suppressMessages({
	    retval =   retval |>
	            dplyr::full_join(object@nodeSS) |>
	            dplyr::mutate(relativeRisk = .data$observed/.data$expected, excessCases = .data$observed - .data$expected) |>
	            dplyr::select(.data$node,.data$observed,.data$expected,.data$relativeRisk, .data$excessCases,.data$LRT,.data$pvalue) 
	  })
	  }
	  return(retval)
	  })




