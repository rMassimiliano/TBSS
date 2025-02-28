#' Bernoulli TBSS 
#' @description This function implements a variety of  Bernoulli TBSS, including stratified analysis that can be used with matched data with variable matching ratios. 

#' @param case a \code{matrix} or \code{data.frame} with two columns (leavesID, cases). Here cases indicate the number of events in the treated group.
#' @param control a \code{matrix} or \code{data.frame} with two columns (leavesID, controls). Here control indicates the number of events in the untreated group.
#' @param tree a \code{matrix} or \code{data.frame} with two columns (nodeID, parentID) that specifies the hierarchical structure in an edge list format.
#' @param p a number in (0,1) indicating the probability of being associated with the treated group under the null hypothesis that treated and untreated have the same probability of experiencing an event. When \code{type = "variableMR"}, this has to be specified in case and control matrices. See details for more information.
#' @param min_events the minimum number of observed events for a node to be tested. This is an integer greater or equal to 2. Default \code{min_events = 2}. Note that by choosing higher \code{min_events}, the power can increase because of the decreased number of hypotheses that are tested. This argument should be fixed before running the analysis to preserve the validity of the p-values. 
#' @param B number of Monte Carlo replicates. The default is \code{B = 999}.
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.
#'@param type determines the type of Bernoulli TBSS to perform. Possible values are \code{type = "unconditional"} for unconditional Bernoulli TBSS, \code{type = "conditional"} for TBSS conditional to the total number of events in the tree, and  \code{type = "variableMR"} for variable matching ratios TBSS (see Details). Default is \code{type = "unconditional"}. 
#'@param nodeToRemove a character vector with user-specified nodes to be removed from the analysis.  Note that sufficient statistics for these nodes are computed, but they are excluded from the testing procedure. By removing nodes from the list of tests, power can increase because of the decreased number of hypotheses that are tested. This list should be specified before running the analysis.  
#'@param seed a number indicating the seed used for the analysis, i.e., the input of \code{set.seed()}. Default \code{seed  = 8277}.
#'@details In TBSS, analyzing matched data with a variable matching ratio requires specifying a different \eqn{p} for each matched set. Typically, we have a length-\eqn{M} vector \eqn{(p_1, \ldots,p_M)} corresponding to \eqn{m:1} matched sets with \eqn{p_m = 1/(1+m)} for \eqn{m=1,\ldots,M.} Note that the same logic applies to any stratified analysis where the \eqn{M} sets are strata, and \eqn{p_m} is the probability of being associated with the untreated group under the null hypothesis. An example of how data are formatted for this analysis is provided in the examples.



#'@export
#' @examples
#' ## Unconditional Bernoulli TBSS examples
#' library(dplyr)
#' library(TBSS)
#'
#' # load tree 
#' data("tree_example", package = 'TBSS')
#'
#' # load data
#' data("bernoulli_TBSS_examples", package = 'TBSS')
#'
#' #treated group
#' data_treated = bernoulli_TBSS_examples$fixed_matching_ratio$case
#' head(data_treated)
#'
#' # control group
#' data_untreated = bernoulli_TBSS_examples$fixed_matching_ratio$control
#' head(data_untreated)
#'
#'matching_prob = bernoulli_TBSS_examples$fixed_matching_ratio$p
#'
#'mod_tbss = bernoulliTBSS(case = data_treated,
#'                         control = data_untreated, 
#'                         tree = tree_example, 
#'                          p = matching_prob )
#' # inspect potential signals
#' summary(mod_tbss) |> filter(pvalue <=0.05)
#'
#'
#' ## Bernoulli TBSS stratified (e.g., for matched data with variable matching ratio) 
#'
#' library(dplyr)
#' library(TBSS)
#'
#' # load tree 
#' data("tree_example", package = 'TBSS')
#'
#' # load data
#' data("bernoulli_TBSS_examples", package = 'TBSS')
#'data_treated = bernoulli_TBSS_examples$variable_matching_ratio$case
#'head(data_treated)
#'
#'data_untreated = bernoulli_TBSS_examples$variable_matching_ratio$control
#'head(data_untreated)
#'
#' mod_tbss = bernoulliTBSS(case = data_treated, control = data_untreated, tree = tree_example, type ='variableMR')
#'
#' # inspect potential signals
#' summary(mod_tbss) |> filter(pvalue <=0.05)

bernoulliTBSS =  function(case,control,tree, p = 1/2, min_events = 2, B = 999, parallel = FALSE, ncpus = NULL, type ='unconditional', nodeToRemove = NULL, seed  = 8277)
{
set.seed(seed)
## stop if the parameter min_events is smaller than 2 or negative 
if(min_events < 2) stop("min_events should be a number greater than 2.")
 ## create data 
 if(type == 'unconditional' | type =='conditional')  
 {
  names(case) = c("leaf", "case")
  names(control) = c("leaf", "control")
  data =  case |> dplyr::full_join(control, by = 'leaf') |> as.data.frame()
  ## NAs are 0s.
  data = data |> dplyr::mutate_at(2:3,\(x) coalesce(x,0L))
 }
 else if(type == 'variableMR')
 {
  names(case) = c("leaf", "case","p")
  names(control) = c("leaf", "control","p")
data =  case |> dplyr::full_join(control, by = c('leaf','p')) |> select('leaf', 'case', 'control','p') |>as.data.frame() 
## check if there are ps that are Nas
 # data$p.x[is.na(data$p.x)] = data$p.y[is.na(data$p.x)]
 #data$p.y[is.na(data$p.y)] = data$p.x[is.na(data$p.y)]
  ## NAs are 0s.

  data = data |> dplyr::mutate_at(2:3,\(x) coalesce(x,0L))
  # data |> dplyr::select(-p.y)
  # names(data)[which(names(data) =='p.x')] = 'p'
 } 
 else stop(sprintf("invalid type %s \n", type))

 tree = as.data.frame(tree)
 if(parallel & (length(ncpus) == 0))
 {
  ncpus = detectCores() - 1
 }

 ## create appropriate TS object
 if(type == 'unconditional')
 {
  myTS = new("unconditionalBernoulliTS", data =data, tree = tree, B = B, p = p)
 }
 else if(type == 'conditional')
 {
   stop("type = 'conditional' is not implemented yet")
 }
 else if(type == 'variableMR')
 {
  myTS = new("unconditionalBernoulliTSVariableRR", data =data, tree = tree, B = B)
 }

  ## TBSS analysis
  myTS = mapNodesLeaves(myTS)
  myTS = computeSS(myTS)

  ## Remove nodes with a number of events lower than min_events 
   myTS@nodesToTest = myTS@nodeSS$node[(myTS@nodeSS$case + myTS@nodeSS$control) >= min_events]
  ## Remove user specified node
  if(!is.null(nodeToRemove))
  {
   myTS@nodesToTest = setdiff(myTS@nodeSS$node,nodeToRemove) 
  }

  myTS = computeLRT(myTS)
  myTS = monteCarlo(myTS, parallel, ncpus)
  return(myTS)
}

