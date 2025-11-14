#############################################
#############################################
### Unconditional Possion Tree Scan
#############################################
#############################################
setClass(Class = "unconditionalPoissonTS",  contains = 'TS')


setMethod("computeSS", signature(object = "unconditionalPoissonTS"), function(object)
{
  object@nodeSS = data.frame(node = names(object@mapNodesLeaves), observed = 0L, expected = 0L)

  for(r in 1:length(object@mapNodesLeaves))
  {

   #ind = which(object@data$leave %in% object@mapNodesLeaves[[r]])
   ind = which(object@leaves %in% object@mapNodesLeaves[[r]])
   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   object@nodeSS$observed[r] = sum(object@data$observed[ind])
   object@nodeSS$expected[r] = sum(object@data$expected[ind])
  }
  return(object)
}
)


setMethod("computeLRT", signature(object = "unconditionalPoissonTS"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## compute LLR for all the nodes
  LRT = numeric(length(object@nodesToTest))
  for(r in 1:length(object@nodesToTest))
  {
   j = which(object@nodeSS$node == object@nodesToTest[r])
   LRT[r] = uncon_poisson_llrt(object@nodeSS$observed[j],object@nodeSS$expected[j])
  }
   object@LRT = LRT
   names(object@LRT) = object@nodesToTest
   return(object)
 })

#' generate data under the null
setMethod("H0_gen", signature(object = "unconditionalPoissonTS"),
function (object)
{
    ## generate the leaves data under H0
    ## we generate data from the leaves even if these are not tested	  
     leaves_list = object@leaves   

     indx = which(object@nodeSS$node %in% leaves_list)
     object@nodeSS$observed[indx] = rpois(length(indx), object@nodeSS$expected[indx])
  
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
	  

setMethod(f = "summary",
          signature = "unconditionalPoissonTS",
          definition = function(object){

          ranking   = order(object@LRT,decreasing = TRUE)
          boot_maxs = apply(object@LRT_H0,2,max)
          pvalues   = sapply(object@LRT[ranking], \(x) (sum(boot_maxs >= x) +1)/(object@B +1))


	  retval = data.frame(node = names(object@LRT)[ranking],
		                     LRT = object@LRT[ranking],
		                     pvalue =  pvalues)
          rownames(retval) = NULL
	  suppressMessages({
	retval =   retval |>
	           full_join(object@nodeSS) |>
	           mutate(relativeRisk = observed/expected, excessCases = observed - expected) |>
	           select(node,observed,expected,relativeRisk, excessCases,LRT,pvalue) 
	  })
	  return(retval)
	  })




setMethod("H1_gen", signature(object = "unconditionalPoissonTS"),
  function (object, alternativeParameter, alternativeLeaves)
  {
    ## generate the leaves data under H0
    ## we generate data from the leaves even if these are not tested	  
     #leaves_list =  object@leaves[which(object@leaves %in% object@nodesToTest)]
     leaves_list = object@leaves   

     ## Null leaves
     indx = which(object@nodeSS$node %in% leaves_list)
     object@nodeSS$observed[indx] = rpois(length(indx), object@nodeSS$expected[indx])

     ## Alternative leaves
     indx = which(object@nodeSS$node %in% alternativeLeaves)
     object@nodeSS$observed[indx] = rpois(length(indx), object@nodeSS$expected[indx]*alternativeParameter)
  
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
