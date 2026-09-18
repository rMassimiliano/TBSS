setClass(Class = "robustSurvivalTS",  contains = 'TS')


setMethod("computeSS", signature(object = "robustSurvivalTS"), function(object, parallel = FALSE, ncpus = NULL)
{
  ## data frame that includes Sufficient Statistics For cases and control
  object@nodeSS = data.frame(node = names(object@mapNodesLeaves),
                    sample_size_case = 0L,
                    events_case = 0L,
                    total_person_time_case = 0L,
                    total_person_time_SOS_case = 0L,
                    uncensored_person_time_case = 0L,
                    sample_size_control = 0L,
                    events_control = 0L,
                    total_person_time_control = 0L,
                    total_person_time_SOS_control = 0L,
                    uncensored_person_time_control = 0L)

  for(r in 1:length(object@mapNodesLeaves))
  {
   ind = which(object@leaves %in% object@mapNodesLeaves[[r]])

   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   
   object@nodeSS$sample_size_case[r] = sum(object@data$sample_size_case[ind])
   object@nodeSS$events_case[r]      = sum(object@data$events_case[ind])
   object@nodeSS$total_person_time_case[r] = sum(object@data$total_person_time_case[ind])
   object@nodeSS$total_person_time_SOS_case[r] = sum(object@data$total_person_time_SOS_case[ind])
   object@nodeSS$uncensored_person_time_case[r] = sum(object@data$uncensored_person_time_case[ind])
   object@nodeSS$sample_size_control[r] = sum(object@data$sample_size_control[ind])
   object@nodeSS$events_control[r] = sum(object@data$events_control[ind])
   object@nodeSS$total_person_time_control[r] = sum(object@data$total_person_time_control[ind])
   object@nodeSS$total_person_time_SOS_control[r] = sum(object@data$total_person_time_SOS_control[ind])
   object@nodeSS$uncensored_person_time_control[r] = sum(object@data$uncensored_person_time_control[ind])
  }
  return(object)
}
)


setMethod("computeLRT", signature(object = "robustSurvivalTS"), function(object) {
  if(length(object@mapNodesLeaves) == 0) object = mapNodesLeaves(object)
  if(length(object@nodeSS) == 0) object = computeSS(object)
  if(length(object@nodesToTest) == 0) object@nodesToTest = names(object@mapNodesLeaves)

  ## compute LLR for all the nodes
  LRT = numeric(length(object@nodesToTest))
  for(r in 1:length(object@nodesToTest))
  {
   j = which(object@nodeSS$node == object@nodesToTest[r])

  # LRT[r] = abs((log(object@nodeSS$events_case[j]) - log(object@nodeSS$total_person_time_case[j]) - (log(object@nodeSS$events_control[j]) - log(object@nodeSS$total_person_time_control[j])))/sqrt(1/object@nodeSS$events_case[j] + 1/object@nodeSS$events_control[j]))
  lambdaH0 = (object@nodeSS$events_case[j] + object@nodeSS$events_control[j])/
(object@nodeSS$total_person_time_case[j] + object@nodeSS$total_person_time_control[j])

LRT[r] = 2*(object@nodeSS$events_case[j]*(log(object@nodeSS$events_case[j]) - log(object@nodeSS$total_person_time_case[j])) - object@nodeSS$events_case[j] + 
(object@nodeSS$events_control[j]*(log(object@nodeSS$events_control[j]) - log(object@nodeSS$total_person_time_control[j])) - object@nodeSS$events_control[j]) -(log(lambdaH0)* (object@nodeSS$events_case[j] + object@nodeSS$events_control[j])  - (object@nodeSS$events_case[j] + object@nodeSS$events_control[j])))



  }
   object@LRT = LRT
   names(object@LRT) = object@nodesToTest
   return(object)
 })


setMethod("H0_gen", signature(object = "robustSurvivalTS"),
function (object)
{
    ## generate the leaves data under H0
    ## we generate data from the leaves even if these are not tested	  

    ## compute necessary quantities
     leaves_list = object@leaves   
     indx = which(object@nodeSS$node %in% leaves_list)

    hat_lambda = cbind(object@nodeSS$events_control[indx]/object@nodeSS$total_person_time_control[indx],
     object@nodeSS$events_case[indx]/object@nodeSS$total_person_time_case[indx])
    
    hat_beta = cbind(-log(hat_lambda[,1]), -log(hat_lambda[,2]) +log(hat_lambda[,1]))
    
    
    minus_H_beta = array(NA, dim = c(length(indx),2,2))
    minus_H_beta[,1,1] = exp(-hat_beta[,1])*object@nodeSS$total_person_time_control[indx] + exp(-rowSums(hat_beta))*object@nodeSS$total_person_time_case[indx]
    
    minus_H_beta[,2,2] = minus_H_beta[,1,2] = minus_H_beta[,2,1] =  exp(-rowSums(hat_beta))*object@nodeSS$total_person_time_case[indx]
    
     J =  apply(minus_H_beta, 1,solve,simplify = FALSE)
    
     A =  array(NA, dim = c(length(indx),2,2))
    
     A[,1,1] = (object@nodeSS$events_case[indx] + object@nodeSS$events_control[indx])+
    object@nodeSS$total_person_time_SOS_control[indx] *exp(-2*hat_beta[,1]) + 
    object@nodeSS$total_person_time_SOS_case[indx] *exp(-2*rowSums(hat_beta)) - 2*exp(-hat_beta[,1]) * object@nodeSS$uncensored_person_time_control[indx] - 
    	2*exp(-rowSums(hat_beta)) * object@nodeSS$uncensored_person_time_case[indx]  
    
    
    
    A[,1,2] = A[,2,1]  = A[,2,2] = object@nodeSS$events_case[indx] + object@nodeSS$total_person_time_SOS_case[indx] * exp(-2*rowSums(hat_beta)) -2*exp(-rowSums(hat_beta))  * object@nodeSS$uncensored_person_time_case[indx]  
    
    var_beta_rob = lapply(1:length(indx), \(i) (J[[i]]%*% A[i,,] %*% J[[i]]))


    ## generate leaf data
    boot_hat_beta = array(0,dim =c(length(indx),2))
    sds = sapply(var_beta_rob, \(x) sqrt(x[1,1]))
    boot_hat_beta[,1] = rnorm(length(indx),hat_beta[,1],sds)
    
    cmeans = (boot_hat_beta[,1] - hat_beta[,1])* sapply(var_beta_rob, \(x) x[1,2]/x[1,1])
    
    csds = sqrt(sapply(var_beta_rob, \(x) x[2,2] - x[1,2]*x[2,1]/x[1,1]))
    boot_hat_beta[,2] = rnorm(length(indx),cmeans,csds)

    object@nodeSS$events_control[indx] = exp(-boot_hat_beta[,1]) * object@nodeSS$total_person_time_control[indx]

    object@nodeSS$events_case[indx] = exp(-rowSums(boot_hat_beta)) * object@nodeSS$total_person_time_case[indx]

# ----------------------
  
     ## update only the SS for the nodes that are used in the LRT
     ##----------
     ind_set = which(!(object@nodeSS$node%in% leaves_list))
     ind_set = ind_set[which(object@nodeSS$node[ind_set] %in% object@nodesToTest)]
     for(r in ind_set)
     {
      leaf = object@mapNodesLeaves[[object@nodeSS$node[r]]]

      ## cases
      object@nodeSS$events_case[r]      = sum(object@nodeSS$events_case[object@nodeSS$node %in% leaf])

      ## controls
      object@nodeSS$events_control[r] = sum(object@nodeSS$events_control[object@nodeSS$node %in% leaf])
    }
    return(object)
   }
)
