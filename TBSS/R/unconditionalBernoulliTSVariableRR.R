#############################################
#############################################
### Unconditional Bernoulli Tree Scan
#############################################
#############################################
setClass(Class = "unconditionalBernoulliTSVariableRR", contains = 'unconditionalBernoulliTS')



setMethod("computeSS", signature(object = "unconditionalBernoulliTSVariableRR"), function(object)
{
  object@nodeSS = data.frame(node = names(object@mapNodesLeaves), case = 0L, control = 0L)

  for(r in 1:length(object@mapNodesLeaves))
  {
   ind = which(object@data$leaf %in% object@mapNodesLeaves[[r]])

   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   object@nodeSS$case[r]     =  sum(object@data$case[ind])
   object@nodeSS$control[r]  = sum(object@data$control[ind])

   ## other statistics
   object@nodeSS$expected[r] = sum((object@data$control[ind] + object@data$case[ind]) * object@data$p[ind])
   object@nodeSS$relativeRisk[r] = sum(object@data$case[ind]/object@data$p[ind])/sum(object@data$control[ind]/(1-object@data$p[ind]))
   object@nodeSS$excessCases[r]  = object@nodeSS$case[r] - object@nodeSS$expected[r]
  }
  return(object)
}
)


setMethod("computeLRT", signature(object = "unconditionalBernoulliTSVariableRR"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## compute LLR for all the nodes
  LRT = numeric(length(object@nodesToTest))
  for(r in 1:length(object@nodesToTest))
  {
    leaves   = object@mapNodesLeaves[[object@nodesToTest[r]]]
    dataNode = object@data[object@data$leaf %in% leaves,]
    case     = tapply(dataNode$case, dataNode$p, sum)
    control  = tapply(dataNode$control, dataNode$p, sum)
    p        = as.numeric(names(case))
    LRT[r]   = sum(uncon_bernoulli_llrt_vectorized(case, control, p))
  }

   names(LRT) = object@nodesToTest
   object@LRT = LRT
   return(object)
 })

#' generate data under the null
setMethod("H0_gen", signature(object = "unconditionalBernoulliTSVariableRR"),
  function (object)
  {
    ## generate the leaves data under H0
    ## we generate data from the leaves even if these are not tested	  
     #leaves_list =  object@leaves[which(object@leaves %in% object@nodesToTest)]

     ## generate leafs value from the data
      totals  = object@data$case + object@data$control
      probs   = object@data$p
      cases   = rbinom(length(totals), size = totals, prob = probs)
      controls = totals - cases
      object@data$case    = cases
      object@data$control = controls

     # update nodes to test
     ind_set = 1:length(object@nodeSS$node)
     ind_set = ind_set[which(object@nodeSS$node[ind_set] %in% object@nodesToTest)]
     for(r in ind_set)
     {
      leaf = object@mapNodesLeaves[[object@nodeSS$node[r]]]
      object@nodeSS$case[r] =  sum(object@data$case[object@data$leaf %in% leaf])
      object@nodeSS$control[r] =  sum(object@data$control[object@data$leaf %in% leaf])
    }
    return(object)
  }
)





setMethod(f = "summary",
          signature = "unconditionalBernoulliTSVariableRR",
          definition = function(object){

          ranking = order(object@LRT,decreasing = TRUE)
          boot_maxs = apply(object@LRT_H0,2,max)
pvalues = sapply(object@LRT[ranking], \(x) (sum(boot_maxs >= x) +1)/(object@B +1))

	  retval = data.frame(node = names(object@LRT)[ranking],
		                     LRT = object@LRT[ranking],
		                     pvalue =  pvalues)

          suppressMessages({
		  retval = retval |> full_join(object@nodeSS) |>
		   select(node,case,control,expected,relativeRisk, excessCases,LRT,pvalue);
		  retval = retval |> arrange(pvalue);
		  rownames(retval) = NULL;
	          })

	  return(retval)
	  })

