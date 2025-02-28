#############################################
#############################################
### Unconditional case case time control Tree Scan
#############################################
#############################################
setClass(Class = "unconditionalCaseCaseTimeControlTS", slots = c(p = 'numeric'), contains = 'TS')

setMethod("computeSS", signature(object = "unconditionalCaseCaseTimeControlTS"), function(object)
{

  object@nodeSS = data.frame(node = names(object@mapNodesLeaves),
                             control_reference = 0,
                             control_hazard = 0,
                             case_reference = 0,
                             case_hazard = 0)


  for(r in 1:length(object@mapNodesLeaves))
  {
   ind = which(object@leaves %in% object@mapNodesLeaves[[r]])
   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   object@nodeSS$case_reference[r] = sum(object@data$case_reference[ind])
   object@nodeSS$case_hazard[r] = sum(object@data$case_hazard[ind])
   object@nodeSS$control_reference[r] = sum(object@data$control_reference[ind])
   object@nodeSS$control_hazard[r] = sum(object@data$control_hazard[ind])
  }
  return(object)
}
)

setMethod("computeLRT", signature(object = "unconditionalCaseCaseTimeControlTS"), function(object){

  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)

  ## compute LLR for all the nodes
  LRT = numeric(length(object@nodesToTest))
  for(r in 1:length(object@nodesToTest))
  {
   j = which(object@nodeSS$node == object@nodesToTest[r])
   LRT[r] = cctc_diff_log_or(object@nodeSS$control_reference[j],
			     object@nodeSS$control_hazard[j],
			     object@nodeSS$case_reference[j], 
			     object@nodeSS$case_hazard[j])

  }
   object@LRT = LRT
   names(object@LRT) = object@nodesToTest
   return(object)
 })


setMethod("H0_gen", signature(object = "unconditionalCaseCaseTimeControlTS"),
function (object)
{
 leaves_list = object@leaves   
 indx = which(object@nodeSS$node %in% leaves_list)

 n_control = object@nodeSS$control_reference[indx] + object@nodeSS$control_hazard[indx]
 n_case = object@nodeSS$case_reference[indx] + object@nodeSS$case_hazard[indx]
  p_null = (object@nodeSS$control_hazard[indx])/(object@nodeSS$control_reference[indx] + object@nodeSS$control_hazard[indx])

   c_case = rbinom(length(indx),n_case, prob = p_null) 	 
   c_control = rbinom(length(indx),n_control, prob = p_null) 	 

   object@nodeSS$control_reference[indx] = n_control - c_control
   object@nodeSS$control_hazard[indx]    = c_control
   object@nodeSS$case_reference[indx]    = n_case - c_case
   object@nodeSS$case_hazard[indx]       = c_case

   ind_set = which(!(object@nodeSS$node%in% leaves_list))
   ind_set = ind_set[which(object@nodeSS$node[ind_set] %in% object@nodesToTest)]
    for(r in ind_set)
    {
      leaf = object@mapNodesLeaves[[object@nodeSS$node[r]]]
      ind = which(object@data$leaf %in% leaf)

      object@nodeSS$case_reference[r] = sum(object@data$case_reference[ind])
      object@nodeSS$case_hazard[r] = sum(object@data$case_hazard[ind])
      object@nodeSS$control_reference[r] = sum(object@data$control_reference[ind])
      object@nodeSS$control_hazard[r] = sum(object@data$control_hazard[ind])
    }

   return(object)
}
)


setMethod(f = "summary",
          signature = "unconditionalCaseCaseTimeControlTS",
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
	           select(node, case_reference, control_reference, case_hazard, control_hazard,LRT,pvalue) })
	  return(retval)
	  })

