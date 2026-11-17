#' Robust survival TBSS 
#' @description This function implements a robustified version of the  exponential survival (or poisson survival) TBSS model that uses exponentially distributed time to event for the leaf level events as a working model rather than the correctly specified one. 
#' @param data a \code{matrix} or \code{data.frame} with 12 columns leavesID, sample_size_case, events_case, total_person_time_case, total_person_time_SOS_case, uncensored_person_time_case,sample_size_control, events_control, total_person_time_control, total_person_time_SOS_control, uncensored_person_time_control), where SOS = Sum of Squares
#' @param tree a \code{matrix} or \code{data.frame} with two columns (nodeID, parentID) that specifies the hierarchical structure in an edge list format.
#' @param min_events the minimum number of observed events for a node to be tested. This is an integer greater or equal to 1. Default \code{min_events = 1}. Note that by choosing higher \code{min_events}, the power can increase because of the decreased number of hypotheses that are tested. This argument should be fixed before running the analysis to preserve the validity of the p-values.  
#' @param B number of Monte Carlo replicates. The default is \code{B = 999}.
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.
#'@param direction can be either \code{'positive'}, \code{"negative"} or \code{"all"}. By default (\code{direction = 'positive'}), the LRT of the nodes where the observed number of cases (or rate) is not larger than the expected one is set to 0. The restriction is applied in the same way to the observed data and to every Monte Carlo replicate. Setting \code{direction = "all"} tests all nodes, while \code{direction = "negative"} only considers nodes where the observed number of cases (or rate) is not larger than the expected.
#'@param nodeToRemove a character vector with user-specified nodes to be removed from the analysis.  Note that sufficient statistics for these nodes are computed, but they are excluded from the testing procedure. By removing nodes from the list of tests, power can increase because of the decreased number of hypotheses that are tested. This list should be specified before running the analysis.  
#'@param seed a number indicating the seed used for the analysis, i.e., the input of \code{set.seed()}. Default \code{seed  = 8277}. If \code{seed = NULL} no seed is used.
#'@references 
#' Russo, M., Hahn, G., Huybrechts, K. F., & Wang, S. V. (2026). A tree-based scan statistic for database studies with time-to-event outcomes. arXiv preprint arXiv:2609.38582.
#'@export
#'@examples 
#' ## Unconditional Poisson TBSS examples
#' library(TBSS)
#'
#' # load tree 
#' data("tree_example", package = 'TBSS')
#' # load data
#' data("survival_poisson_example", package = 'TBSS')
#' dat =  survival_poisson_example 
#' head(dat)
#'
#'
#'
#'mod_tbss = robustSurvivalTBSS(data = dat, tree = tree_example, direction = 'all')
#' # inspect potential signals
#' summary(mod_tbss) |> dplyr::filter(pvalue <=0.05)
robustSurvivalTBSS =  function(data,tree, min_events = 1, B = 999, parallel = FALSE, ncpus = NULL, direction = 'positive', nodeToRemove = NULL, seed  = 8277)
{
 if(!is.null(seed)) set.seed(seed)
 ## stop if the parameter min_events is smaller than 1 or negative 
 if(min_events < 1) stop("min_events should be a number greater than 1.")
 if(!(direction %in% c("positive","all","negative"))) stop("direction should be 'positive', 'all', or 'negative'. Refer to the help for details.")

 names(data) = c("leaf", "sample_size_case", "events_case", "total_person_time_case", "total_person_time_SOS_case", "uncensored_person_time_case", "sample_size_control", "events_control", "total_person_time_control", "total_person_time_SOS_control",  "uncensored_person_time_control")

 tree = as.data.frame(tree)
 if(parallel & (length(ncpus) == 0))
 {
  ncpus = parallel::detectCores() - 1
 }

 ## Create appropriate TS object
   myTS = methods::new("robustSurvivalTS", tree = as.data.frame(tree), data = as.data.frame(data), B = B)
   myTS@direction = direction
  ## TBSS analysis
  myTS = mapNodesLeaves(myTS)
  myTS = computeSS(myTS)

  ## Remove nodes with a number of events lower than min_events 
  myTS@nodesToTest = myTS@nodeSS$node[(myTS@nodeSS$events_case >= min_events) & (myTS@nodeSS$events_control>0)] 

  ## Remove user specified node
  if(!is.null(nodeToRemove))
  {
   myTS@nodesToTest = setdiff(myTS@nodesToTest, nodeToRemove) 
  }
  if(length(myTS@nodesToTest) == 0) stop("No outcome node to be tested after removing nodes from the list for min_events conditions and explicitly with nodeToRemove.")

  ## compute node-specific discrepancies and their distribution under the null
  myTS = computeLRT(myTS)
  myTS = monteCarlo(myTS, parallel, ncpus)
  return(myTS)
}

