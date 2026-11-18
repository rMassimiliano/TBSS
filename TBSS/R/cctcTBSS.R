#' This function implements TreeScan for a case-case-time control design
#' @param `case` a matrix or data.frame with three columns (leavesID, cases in the reference window, cases in the hazard window)
#' @param `control` a matrix or data.frame with three columns (leavesID, cases in the reference window, cases in the hazard window)
#' @param `tree` a matrix or data.frame with two columns (nodeID, parentNodeID) 
#' @param `B` number of bootstrap replicates. Default is `999`.
#' @param `parallel` if TRUE parallel computing is used to compute the bootstrap p-values.
#' @param `ncpus` number of cpus used if `parallel = TRUE`. If no value is provided set to  `detectCores() - 1`.
#'@param nodeToRemove a character vector with user-specified nodes to be removed from the analysis.  Note that sufficient statistics for these nodes are computed, but they are excluded from the testing procedure. By removing nodes from the list of tests, power can increase because of the decreased number of hypotheses that are tested. This list should be specified before running the analysis. 
#'@param seed a number indicating the seed used for the analysis, i.e., the input of \code{set.seed()}. Default \code{seed  = 8277}.If \code{seed = NULL} no seed is used.
#' @return An S4 object of class \code{unconditionalCaseCaseTimeControlTS}. The slot \code{LRT} contains the log-likelihood ratio of each tested node (named by node). The tested nodes are also listed in slot \code{nodesToTest}). The slot \code{LRT_H0} is a matrix with one row per tested node and one column per Monte Carlo replicate, and the slot \code{nodeSS} contains the sufficient statistics of every node of the tree. Use \code{summary()} for a table of the results with the Monte Carlo p-values.
#' @export
#'@examples 
#' ## Case-Case-Time-Control example
#' library(TBSS)
#'
#' # load tree 
#' data("tree_example", package = 'TBSS')
#'
#' # load data
#' data("cctc_TBSS_example", package = 'TBSS')
#'
#'
#'mod_tbss = cctcTBSS(case =cctc_TBSS_example$case, control = cctc_TBSS_example$control,  tree = tree_example)
#' summary(mod_tbss) 
cctcTBSS =  function(case,control,tree,B = 999, parallel = FALSE, ncpus = NULL, nodeToRemove = NULL, seed  = 8277)
{
 if(!is.null(seed)) set.seed(seed)
 ## create data 
 names(case) = c("leaf", "case_reference", "case_hazard")
 names(control) = c("leaf", "control_reference", "control_hazard")
 data = control |> dplyr::full_join(case, by = 'leaf') |> as.data.frame()
 ## replace NA with 0s (for leaves with only cases or controls)
 data = data |> dplyr::mutate_at(2:5,\(x) dplyr::coalesce(x,0L))

 tree = as.data.frame(tree)
 if(parallel & (length(ncpus) == 0))
 {
  ncpus = parallel::detectCores() - 1
 }

 ## create appropriate TS object
  myTS = methods::new("unconditionalCaseCaseTimeControlTS", data =data, tree = tree, B =B )
  myTS = mapNodesLeaves(myTS)
  myTS = computeSS(myTS)
  myTS@nodesToTest = myTS@nodeSS$node
  ## Remove user specified node
  if(!is.null(nodeToRemove))
  {
   myTS@nodesToTest = setdiff(myTS@nodesToTest, nodeToRemove) 
  }
  if(length(myTS@nodesToTest) == 0) stop("No outcome node to be tested after removing nodes from the list for min_events conditions and explicitly with nodeToRemove.")
  myTS = computeLRT(myTS)
  myTS = monteCarlo(myTS, parallel, ncpus)

  return(myTS)
}

