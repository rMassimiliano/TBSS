#' Poisson  log-likelihood ratio statistic
#' compute the log-likelihood ratio statistic for the difference in mean of two Poisson distributions
#' @param observed number of events for the first Poisson
#' @param expected number of events for the second Poisson
uncon_poisson_llrt =function (observed,expected)
{
   if(observed==0){return(expected)}
   res = expected - observed + observed*(log(observed) - log(expected)  )
   return(res)
}

#' Poisson  log-likelihood ratio statistic
#' compute the log-likelihood ratio statistic for  for the difference in mean of two Poisson distributions with different offsets (person time)
#' @param observed number of events for the first Poisson
#' @param observed_pt  time reference (or person time) for the first Poisson
#' @param expected number of events for the second Poisson
#' @param expected_pt time reference (or person time) for the second Poisson
poisson_llrt =function (observed, observed_pt, expected, expected_pt)
{
   tot = (expected + observed)
   l0 = tot/(observed_pt + expected_pt)

   res =  ifelse(observed>0,-observed + observed*log(observed),0) +
   ifelse(expected>0, -expected + expected*log(expected),0) + 
   tot - tot*log(l0) - observed*log(observed_pt) - expected * log(expected_pt)
   return(res)
}
