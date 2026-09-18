setClass(Class = "poissonSurvivalTS",  contains = 'TS')


setMethod("computeSS", signature(object = "poissonSurvivalTS"), function(object, parallel = FALSE, ncpus = NULL)
{
  ## data frame that includes Sufficient Statistics For cases and control
  object@nodeSS = data.frame(node = names(object@mapNodesLeaves),
                    events_case = 0L,
                    total_person_time_case = 0L,
                    events_control = 0L,
                    total_person_time_control = 0L)

  for(r in 1:length(object@mapNodesLeaves))
  {
   ind = which(object@leaves %in% object@mapNodesLeaves[[r]])

   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   
   object@nodeSS$events_case[r]      = sum(object@data$events_case[ind])
   object@nodeSS$total_person_time_case[r] = sum(object@data$total_person_time_case[ind])
   object@nodeSS$events_control[r] = sum(object@data$events_control[ind])
   object@nodeSS$total_person_time_control[r] = sum(object@data$total_person_time_control[ind])
  }
  return(object)
}
)


setMethod("computeLRT", signature(object = "poissonSurvivalTS"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## compute LLR for all the nodes
  LRT = numeric(length(object@nodesToTest))
  for(r in 1:length(object@nodesToTest))
  {
   j = which(object@nodeSS$node == object@nodesToTest[r])
   y = c(object@nodeSS$events_case[j], object@nodeSS$events_control[j])
   o = c(object@nodeSS$total_person_time_case[j],  object@nodeSS$total_person_time_control[j])
   LRT[r] = 2*(sum(dpois(y,y, log = T)) - sum( dpois(y, sum(y)/sum(o) *o, log = T)))
  }
   object@LRT = LRT
   names(object@LRT) = object@nodesToTest
   return(object)
 })


setMethod("H0_gen", signature(object = "poissonSurvivalTS"),
function (object)
{
    ## generate the leaves data under H0
    ## we generate data from the leaves even if these are not tested	  
     leaves_list = object@leaves   
     indx = which(object@nodeSS$node %in% leaves_list)

       ## leaf MLEs under the null 
       lH0 = with(object@nodeSS,
(events_case[indx] + events_control[indx])/
     (total_person_time_case[indx] + total_person_time_control[indx]))


     ## generate total events for the leaves in the two gropus
 
     ## generate the total person time from the asymptotic conditional distribution
     ## cases 
     ##----------

       r_case   =   with(object@nodeSS, rpois(length(indx), lH0*  total_person_time_case[indx]))
       object@nodeSS$events_case[indx]  = r_case
     ## controls
     ##----------
       r_control   =   with(object@nodeSS, rpois(length(indx), lH0 * total_person_time_control[indx]))
       object@nodeSS$events_control[indx]  = r_control

  
     ## update only the SS for the nodes that are used in the LRT
     ##----------
     ind_set = which(!(object@nodeSS$node%in% leaves_list))
     ind_set = ind_set[which(object@nodeSS$node[ind_set] %in% object@nodesToTest)]
     for(r in ind_set)
     {
      leaf = object@mapNodesLeaves[[object@nodeSS$node[r]]]

      ## cases
      object@nodeSS$events_case[r]      = sum(object@nodeSS$events_case[object@nodeSS$node %in% leaf])
      object@nodeSS$total_person_time_case[r] = sum(object@nodeSS$total_person_time_case[object@nodeSS$node %in% leaf])

      ## controls
      object@nodeSS$events_control[r] = sum(object@nodeSS$events_control[object@nodeSS$node %in% leaf])
      object@nodeSS$total_person_time_control[r] = sum(object@nodeSS$total_person_time_control[object@nodeSS$node %in% leaf])

    }
    return(object)
   }
)
