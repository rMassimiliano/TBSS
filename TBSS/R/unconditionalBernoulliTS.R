#############################################
#############################################
### Unconditional Bernoulli Tree Scan
#############################################
#############################################


setMethod("computeSS", signature(object = "unconditionalBernoulliTS"), function(object)
{
  object@nodeSS = data.frame(node = names(object@mapNodesLeaves), case = 0L, control = 0L)

  for(r in 1:length(object@mapNodesLeaves))
  {
   #ind = which(object@data$leaf %in% object@mapNodesLeaves[[r]])

   ind = which(object@leaves %in% object@mapNodesLeaves[[r]])
   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   object@nodeSS$case[r]     =  sum(object@data$case[ind])
   object@nodeSS$control[r]  = sum(object@data$control[ind])
   
  }
  return(object)
}
)

setMethod("computeLRT", signature(object = "unconditionalBernoulliTS"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## compute LLR for all the nodes
  LRT = numeric(length(object@nodesToTest))
  for(r in 1:length(object@nodesToTest))
  {
   j = which(object@nodeSS$node == object@nodesToTest[r])
   LRT[r] =  uncon_bernoulli_llrt(object@nodeSS$case[j],object@nodeSS$control[j],object@p)
  }
   names(LRT) = object@nodesToTest
   object@LRT = LRT
   return(object)
 })


#' generate a new tree under the global null H0
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
     cases  = rbinom(length(indx), size = totals, prob = object@p)
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
		  retval = retval |> full_join(object@nodeSS);
                  retval = retval |>
	           mutate(expected = (control +case)*object@p) |>
		   mutate(relativeRisk =  case/expected,
			  excessCases = case - expected) |>
		   select(node,case,control, expected, relativeRisk, excessCases, LRT,pvalue);
		  retval = retval |> arrange(pvalue);
		  rownames(retval) = NULL;
	          })

	  return(retval)
	  })

