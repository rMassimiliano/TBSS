#' @importFrom rlang .data
#' @importFrom survival Surv coxph
setClass(Class = "coxSurvivalTS",  contains = 'permutationPoissonTS')

#' Sufficient Statistics
#' Compute sufficient statistics for the permutationPoissonTS 
#' @param object an object of class TS
#' @rdname permutationPoissonTS
setMethod("computeSS", signature(object = "coxSurvivalTS"), function(object)
{
  object@nodeSS = data.frame(node = names(object@mapNodesLeaves), observed = 0L, observed_pt = 1L, expected = 0L, expected_pt = 1L, HR = 0L)

  for(r in 1:length(object@mapNodesLeaves))
  {

   object@nodeSS$node[r]     = names(object@mapNodesLeaves)[r]
   object@nodeSS$observed[r] = sum(object@idInNode[[r]]$weight*(object@idInNode[[r]]$exposure == 1))
   object@nodeSS$expected[r] = sum( object@idInNode[[r]]$weight* (object@idInNode[[r]]$exposure == 0))
    ## combine person time for censored and not in node
    combined_pt = rbind(object@personTime[[r]], object@personTimeNotInNode[[r]])

    object@nodeSS$expected_pt[r] = with(dplyr::filter(object@idInNode[[r]],.data$exposure == 0), sum(weight*followup)) + with(dplyr::filter(combined_pt,.data$exposure == 0), sum(weight*followup)) 

    object@nodeSS$observed_pt[r] =  with(dplyr::filter(object@idInNode[[r]],.data$exposure == 1), sum(weight*followup)) + with(dplyr::filter(combined_pt,.data$exposure == 1), sum(weight*followup)) 
   
    ## estimate cox model

	t2 = object@idInNode[[r]] |> dplyr::select(.data$patientID,.data$followup,.data$exposure,.data$weight) |> dplyr::mutate(status=1)  |> dplyr::rename(time =  .data$followup)
  	t1 = object@personTime[[r]] |> dplyr::select(.data$patientID,.data$followup,.data$exposure,.data$weight)|> dplyr::mutate(status = 0)|> dplyr::rename(time =.data$followup)
  	t3 = object@personTimeNotInNode[[r]] |> dplyr::select(.data$patientID,.data$followup,.data$exposure,.data$weight)|> dplyr::mutate(status = 0)|> dplyr::rename(time =.data$followup)
  
  	tmpData  = rbind(t1,t2,t3)  
  
  	if(NROW(tmpData)>2)
  	{
  		cmod = survival::coxph(survival::Surv(time,status)~exposure, weight = tmpData$weight, data = tmpData, ties = "breslow")
		object@nodeSS$HR[r] = exp(cmod$coef[[1]])
  
  	}
	else
	{
		object@nodeSS$HR = NA
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

  d = length(object@nodesToTest)
  LRT = numeric(d)
  coxTBSS = data.frame(node = character(d), coef = numeric(d), LRT = numeric(d))
  i=1
  for(node in object@nodesToTest)
  {
  	coxTBSS$node[i] = node
  	t2 = object@idInNode[[node]] |> dplyr::select(.data$patientID,.data$followup,.data$exposure,.data$weight) |> dplyr::mutate(status=1)  |> dplyr::rename(time = .data$followup)
  	t1 = object@personTime[[node]] |> dplyr::select(.data$patientID, .data$followup, .data$exposure, .data$weight)|> dplyr::mutate(status = 0)|> dplyr::rename(time = .data$followup)
  	t3 = object@personTimeNotInNode[[node]] |> dplyr::select(.data$patientID,.data$followup,.data$exposure,.data$weight)|> dplyr::mutate(status = 0)|> dplyr::rename(time = .data$followup)
  
  
  	tmpData  = rbind(t1,t2,t3)  
  
  	if(NROW(tmpData)>2)
  	{
  		cmod = survival::coxph(survival::Surv(time,status)~exposure, weight = tmpData$weight, data = tmpData, ties = "breslow")
  		LRT[i] = -2 * (cmod$loglik[1] - cmod$loglik[2])
  
  	}
  	i=i+1
  }
   object@LRT = LRT
   names(object@LRT) = object@nodesToTest
   return(object)
 })




#'  Summary method
#' @description provide a  summary of the TBSS analysis results. Include Hazard Ratio (HR), time and relative measures 
#' @param object an object of class TS
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
