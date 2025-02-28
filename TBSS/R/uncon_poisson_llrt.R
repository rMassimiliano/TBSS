#' compute the log-likelihood ratio test for unconditional Poisson analysis 
uncon_poisson_llrt =function (observed,expected)
{
 # res =   dpois(observed,observed, log = TRUE) - dpois(observed,expected, log = TRUE)
   if(observed==0){return(-expected)}
   res = expected - observed + observed*(log(observed) - log(expected)  )
   return(res)
}
