#' compute the log-likelihood ratio test for unconditional binomial analysis 
uncon_bernoulli_llrt =function (case,control, p)
{
  res =  if(case/(case + control) > p & (control >0))
    case*(log(case) - log(p)) - (case+control)*log(case + control) + control*(log(control) - log(1-p)) 
  else if (case/(case + control) > p & (control ==0)) 
	  -case*log(p)
    else 0
  res
}



#' same function but vectorized. 
uncon_bernoulli_llrt_vectorized =function (case,control, p)
{
   res = rep(0,length(case))
   conditions = (case/(case + control) > p) & (control >0)

   Ccase = case[conditions]
   Ccontrol= control[conditions]
   Cp  =  if(length(p)>1)  p[conditions] else rep(p,sum(conditions))

   res[conditions] =  Ccase*(log(Ccase) - log(Cp)) - (Ccase+Ccontrol)*log(Ccase + Ccontrol) + Ccontrol*(log(Ccontrol) - log(1-Cp))
   res
}
