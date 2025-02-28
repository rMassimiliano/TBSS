#' Discrepancy measure for the case case time controll TBSS
#' @description Compute the discrepancy measure used for the case-case-time-control TBSS. The discrepancty is the log difference of the odds ratios for the case-case-time-control analysis. 
#' @param control_reference: number of control in the reference window
#' @param control_hazard: number of control in the hazard window
#' @param case_reference: number of cases in the reference window
#' @param case_hazard: number of cases in the hazard window
cctc_diff_log_or = function(control_reference,
			    control_hazard,
			    case_reference,
			    case_hazard)
{
 	

  LOR_1 = log(case_hazard) - log(case_reference)
  LOR_0 = log(control_hazard) - log(control_reference)
  retval = LOR_1 - LOR_0 
  ## only if positive
  if(is.nan(retval)) retval = 0 ## Inf - Inf !! 
  ## we can use this for unilateral test
  ##res = if(retval >0) retval else 0
  return(retval)
}

