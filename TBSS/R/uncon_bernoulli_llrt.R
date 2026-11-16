
#' log-likelihood ratio for binomial distribution 
#' compute the log-likelihood ratio test statistic for unconditional binomial analysis 
#' @param case number of case (successes)
#' @param control number of control (failures)
#' @param p theoretical probability of the binomial
uncon_bernoulli_llrt =function (case,control, p)
{
  res =  if(case/(case + control) > p & (control >0))
    case*(log(case) - log(p)) - (case+control)*log(case + control) + control*(log(control) - log(1-p)) 
  else if (case/(case + control) > p & (control ==0)) 
	  -case*log(p)
    else 0
  res
}




#' log-likelihood ratio for binomial distribution 
#' compute the log-likelihood ratio test statistic for unconditional binomial analysis (vectorized)
#' @param case number of case (successes)
#' @param control number of control (failures)
#' @param p theoretical probability of the binomial
uncon_bernoulli_llrt_vectorized =function (case,control, p)
{
   res = rep(0,length(case))
   conditions = (case/(case + control) > p) & (control >0)
   conditions0 = (case > 0) & (control == 0)

   Ccase = case[conditions]
   Ccontrol= control[conditions]
   Cp  =  if(length(p)>1)  p[conditions] else rep(p,sum(conditions))

   res[conditions] =  Ccase*(log(Ccase) - log(Cp)) - (Ccase+Ccontrol)*log(Ccase + Ccontrol) + Ccontrol*(log(Ccontrol) - log(1-Cp))
   Cp0 =  if(length(p)>1)  p[conditions0] else rep(p,sum(conditions0))
   res[conditions0] = -case[conditions0]*log(Cp0)
   res
}
