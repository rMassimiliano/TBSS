#' Evaluate power of TBSS for a specific alternative distribution

#' @param object an object of class \code{TS} 
#' @param alternativeParameter either a number or a vector of lenght(alternativeLeaves) specifying the parameters to be used for the alternative model. For Poisson scan statistics this is an increase in relative risk, while for Bernoulli scan statistics is the probability under the alternative hypothesis. 
#' @param alternativeLeaves a vector with the name of the leaves that will be generated under the alternative model. See details
#' @param alpha is the type-I error used for the test. The default is 0.05.
#' @param R is the number of dataset simulated under the alternative and used to estimate the power. 
#' @param parallel if \code{TRUE} parallel computing is used to compute the bootstrap p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided set to  \code{ncpus = detectCores() - 1}.

#' @export
evalPower = function(object,
		     alternativeParameter = 1,
		     alternativeLeaves = NULL,
		     alpha = 0.05,
		     R =1000,
                     parallel = FALSE,
                     ncpus = NULL)
{
 if(length(object@LRT) == 0) object  = computeLRT(object)
 if(length(object@LRT_H0) == 0) object  = monteCarlo(object,parallel,ncpus)

 T_boot = c(apply(object@LRT_H0,2,max), max(object@LRT))
 T_decision = quantile(T_boot, 1-alpha)

 ## function that compute one Monte Carlo replicate for power
 do_MC_pow = function(...)
 {
  boot_object = object
  ## re-generate data under the alternative
  boot_object = H1_gen(object,alternativeParameter, alternativeLeaves)
  statistic = max(TBSS:::computeLRT(boot_object)@LRT)
  return(statistic)
}

if(parallel)
{
 suppressMessages(require(parallel))
 suppressMessages(require(TBSS))
 cl <- makeCluster(ncpus)
 clusterExport(cl, "do_MC_pow", envir = environment())
 clusterExport(cl, "object", envir = environment())
 clusterExport(cl, "alternariveParameter", envir = environment())
 clusterExport(cl, "alternativeLeaves", envir = environment())
 T_alternative = parSapply(cl, 1:R, do_MC_pow)
 stopCluster(cl)
}
else
{
 ## non parallel version
 T_alternative =  sapply(1:R,do_MC_pow)
}
 pow = mean(T_alternative >= T_decision)
 return(pow)
}
