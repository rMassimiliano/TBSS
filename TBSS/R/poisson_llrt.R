
#' compute the log-likelihood ratio test for two Poisson with different offsets (person time)
poisson_llrt =function (observed, observed_pt, expected, expected_pt)
{
 # res =   dpois(observed,observed, log = TRUE) - dpois(observed,expected, log = TRUE)
   tot = (expected + observed)
   l0 = tot/(observed_pt + expected_pt)

   res =  ifelse(observed>0,-observed + observed*log(observed),0) +
   ifelse(expected>0, -expected + expected*log(expected),0) + 
   tot - tot*log(l0) - observed*log(observed_pt) - expected * log(expected_pt)
   return(res)
}
