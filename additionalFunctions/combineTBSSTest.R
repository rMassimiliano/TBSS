#'CombineTBSSTest

#' @description This function combine the results for LRT of a object of class TBSS with LRT from gaussian test. TBSS and Gaussian tests are performed independently but with common multiplicity control.
#' @param tbss an object of class TS
#' @param lab_values a data frame with 5 columns: nodeID, sample size, mean, standard deviations, exposure status
#'@param statistic used on the continuous data. To data only LRT is implemented 
#'@param seed seed for random generation (defulat = 123).

combineTBSSTest = function(tbss, lab_values, statistic = 'LRT', seed = 123) 
{
  set.seed(seed)
  if(statistic != "LRT") stop("Currently the function only support LRT")
  ## set names of lab values
  names(lab_values)[1:5] = c("node","n", "mean", "sd", "exposure")
  
  ## for each node compute the SS & LRT
  SS = lab_values |> 
  group_by(node) |>
  mutate(meanH0 = sum(n*mean)/sum(n))

  tmp = SS |>
  summarize(lrH1 = sum(dnorm(mean,mean,sd/sqrt(n), log = TRUE)),
  lrH0 = sum(dnorm(mean,meanH0,sd/sqrt(n), log = TRUE)),
  LRT = lrH1 - lrH0) 
  
  cLRT  = tmp |> pull(LRT)
  names(cLRT) = tmp$node

  ## generate LRT samples under the null

  bSS = SS
  cLRT_H0 = matrix(NA,length(cLRT),tbss@B)
  
  for(b in seq(tbss@B))
  {
   bSS$mean = rnorm(SS$n,SS$meanH0,SS$sd/sqrt(SS$n))
   cLRT_H0[,b] = 
   bSS |> 
   group_by(node) |>
   mutate(meanH0 = sum(n*mean)/sum(n)) |>
   summarize(lrH1 = sum(dnorm(mean,mean,sd/sqrt(n), log = TRUE)),
   lrH0 = sum(dnorm(mean,meanH0,sd/sqrt(n), log = TRUE)),
   LRT = lrH1 - lrH0) |> pull(LRT)
  }
  ## Combine LRT and Statistics with the Original TBSS

  nLRT = c(tbss@LRT, cLRT)
  nLRT_H0 = rbind(tbss@LRT_H0, cLRT_H0)



  ranking   = order(nLRT,decreasing = TRUE)
  boot_maxs = apply(nLRT_H0,2,max, na.rm = TRUE)
  pvalues   = sapply(nLRT[ranking], \(x) (sum(boot_maxs >= x) +1)/(tbss@B +1))


  retval = data.frame(node = names(nLRT)[ranking],
		                     LRT = nLRT[ranking],
		                     pvalue =  pvalues)
          rownames(retval) = NULL

  finalTable = summary(tbss) |> select(-LRT,-pvalue) |> full_join(retval) |> arrange(pvalue) 
  ## Produce and updated summary tables
  results = list(summary = finalTable, LRT = nLRT, LRT_H0 = nLRT_H0)
  return(results)
}

