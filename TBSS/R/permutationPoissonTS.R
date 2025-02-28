#############################################
#############################################
### Unconditional Possion Tree Scan
#############################################
#############################################
setClass(Class = "permutationPoissonTS",  contains = 'TS', slots = c(idInNode = 'list',personTime = 'list', cohort = 'data.frame', fup = 'logical')) 


setMethod("computeSS", signature(object = "permutationPoissonTS"), function(object)
{
  object@nodeSS = data.frame(node = names(object@mapNodesLeaves), observed = 0L, observed_pt = 1L, expected = 0L, expected_pt = 1L)

  for(r in 1:length(object@mapNodesLeaves))
  {

   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   object@nodeSS$observed[r] = sum(object@idInNode[[r]]$weight*(object@idInNode[[r]]$exposure == 1))
   object@nodeSS$expected[r] = sum( object@idInNode[[r]]$weight* (object@idInNode[[r]]$exposure == 0))
   if(object@fup)
   {
        object@nodeSS$expected_pt[r] = with(dplyr::filter(object@idInNode[[r]],exposure == 0), sum(weight*followup)) + with(dplyr::filter(object@personTime[[r]],exposure == 0), sum(weight*followup)) 


        object@nodeSS$observed_pt[r] =  with(dplyr::filter(object@idInNode[[r]],exposure == 1), sum(weight*followup)) + with(dplyr::filter(object@personTime[[r]],exposure == 1), sum(weight*followup)) 
   }

  }

  
  return(object)
}
)


setMethod("computeLRT", signature(object = "permutationPoissonTS"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## compute LLR for all the nodes
  LRT = numeric(length(object@nodesToTest))
  for(r in 1:length(object@nodesToTest))
  {
   j = which(object@nodeSS$node == object@nodesToTest[r])
   LRT[r] = with(object@nodeSS[j,], poisson_llrt(observed, observed_pt, expected, expected_pt))
  }
   object@LRT = LRT
   names(object@LRT) = object@nodesToTest
   return(object)
 })

#' generate data under the null via permutation of the lables
setMethod("H0_gen", signature(object = "permutationPoissonTS"),
function (object)
{
    ## generate the leaves data under H0
    ## we generate data from the leaves even if these are not tested	  
  #permuted_cohort = object@cohort |> select(strata,exposure,weight) |>  group_by(strata,exposure) |> slice_sample(prop = 1)

	## it is ordered by strata
  permuted_cohort = object@cohort |> dplyr::select(patientID,strata) |>  dplyr::group_by(strata) |> dplyr::slice_sample(prop = 1)

 ## cohort is ordered by strata so this change allocation but not structure since the PatientId are permuted withing the strata
 permuted_cohort$exposure = object@cohort$exposure
 permuted_cohort$weight = object@cohort$weight

     ## update node that have to be tested
     for(j in 1:length(object@nodesToTest))
     {
        r = which(names(object@idInNode) == object@nodesToTest[j])
     ## change exposure and weights for the SS order is given by nodes in idNode
     tmp_cohort = permuted_cohort[permuted_cohort$patientID %in% object@idInNode[[r]]$patientID,]

        object@idInNode[[r]]$exposure  = tmp_cohort$exposure
        object@idInNode[[r]]$weight =   tmp_cohort$weight
      if(object@fup)
      {
        tmp_cohort = permuted_cohort[permuted_cohort$patientID %in% object@personTime[[r]]$patientID,]
      	object@personTime[[r]]$exposure =tmp_cohort$exposure
      	object@personTime[[r]]$weight = tmp_cohort$weight
      object@personTime[[r]]$wft = with(object@personTime[[r]], weight*followup)
      }
     }
     ## compute SS 
     object = computeSS(object)
    return(object)
   }
)
	  

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
	           full_join(object@nodeSS) |>
	           mutate(relativeRisk = observed/expected, excessCases = observed - expected) |>
	           select(node,observed, expected,relativeRisk, excessCases,LRT,pvalue) 
	  })
	  }
	  else
	  {
	   suppressMessages({
	    retval =   retval |>
	           full_join(object@nodeSS) |>
	           mutate(relativeRisk = observed/expected, excessCases = observed - expected) |>
	           select(node,observed,expected,relativeRisk, excessCases,LRT,pvalue) 
	  })
	  }
	  return(retval)
	  })




