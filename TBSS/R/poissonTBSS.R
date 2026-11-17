#' Poisson TBSS 
#' @description This function implements unconditional and conditional Poisson TBSS analysis.
#' @param data a \code{matrix} or \code{data.frame} with three columns (leavesID, observed, expected). 
#' @param tree a \code{matrix} or \code{data.frame} with two columns (nodeID, parentID) that specifies the hierarchical structure in an edge list format.
#' @param min_events the minimum number of observed events for a node to be tested. This is an integer greater or equal to 2. Default \code{min_events = 2}. Note that by choosing higher \code{min_events}, the power can increase because of the decreased number of hypotheses that are tested. This argument should be fixed before running the analysis to preserve the validity of the p-values.  
#' @param B number of Monte Carlo replicates. The default is \code{B = 999}.
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.
#'@param type determines the type of Poisson TBSS to perform. Possible values are \code{type = "unconditional"} for unconditional Poisson TBSS, and \code{type = "conditional"} for Poisson TBSS conditional to the total number of events in the tree.  Default is \code{type = "unconditional"}. 
#'@param direction can be either \code{'positive'}, \code{"negative"} or \code{"all"}. By default (\code{direction = 'positive'}), the LRT of the nodes where the observed number of cases (or rate) is not larger than the expected one is set to 0. The restriction is applied in the same way to the observed data and to every Monte Carlo replicate. Setting \code{direction = "all"} tests all nodes, while \code{direction = "negative"} only considers nodes where the observed number of cases (or rate) is not larger than the expected.
#'@param nodeToRemove a character vector with user-specified nodes to be removed from the analysis.  Note that sufficient statistics for these nodes are computed, but they are excluded from the testing procedure. By removing nodes from the list of tests, power can increase because of the decreased number of hypotheses that are tested. This list should be specified before running the analysis. 
#'@param seed a number indicating the seed used for the analysis, i.e., the input of \code{set.seed()}. Default \code{seed  = 8277}.If \code{seed = NULL} no seed is used.
#'@references
#' Kulldorff, M., Fang, Z. and Walsh, S.J., 2003. A tree-based scan statistic for database disease surveillance. Biometrics, 59(2), pp.323-331.
#'
#' Russo, M. and Wang, S.V., 2024. An open‐source implementation of tree‐based scan statistics. Pharmacoepidemiology and drug safety, 33(3), p.e5765.
#'@export
#'@examples 
#' ## Unconditional Poisson TBSS examples
#' library(TBSS)
#'
#' # load tree 
#' data("tree_example", package = 'TBSS')
#'
#' # load data
#' data("poisson_TBSS_example", package = 'TBSS')
#'
#' dat = poisson_TBSS_example
#' head(dat)
#'
#'
#'
#'mod_tbss = poissonTBSS(data = dat,
#'                         tree = tree_example)
#' # inspect potential signals
#' summary(mod_tbss) |> dplyr::filter(pvalue <=0.05)
#'
# # consider also node with expected < observed
#'mod_tbss = poissonTBSS(data = dat,
#'                       tree = tree_example,
#'                       direction = 'all')
#' # inspect potential signals
#' summary(mod_tbss) |> dplyr::filter(pvalue <=0.05)
#'
#' # conditional analysis
#' mod_tbss_cond = poissonTBSS(data = dat,
#'                         tree = tree_example,
#'                         type = 'conditional'  )
#' # inspect potential signals
#' summary(mod_tbss_cond) |> dplyr::filter(pvalue <=0.05)
poissonTBSS =  function(data,tree, min_events = 2, B = 999, parallel = FALSE, ncpus = NULL, type ='unconditional', direction = 'positive', nodeToRemove = NULL, seed  = 8277)
{
if(!is.null(seed)) set.seed(seed)
## stop if the parameter min_events is smaller than 2
if(min_events < 2) stop("min_events should be a number greater than or equal to 2.")
if(!(direction %in% c("positive","all","negative"))) stop("direction should be 'positive', 'all', or 'negative'. Refer to the help for details.")

 ## create data 
 if(type == 'unconditional' | type =='conditional')  
 {
  names(data) = c("leaf",  "observed", "expected")
 }
 else stop(sprintf("invalid type %s \n", type))

 tree = as.data.frame(tree)
 if(parallel & (length(ncpus) == 0))
 {
  ncpus = parallel::detectCores() - 1
 }

 ## Create appropriate TS object
 if(type == 'unconditional')
 {
   myTS = methods::new("unconditionalPoissonTS", tree = as.data.frame(tree), data = as.data.frame(data), B = B)
 }
 else if(type == 'conditional')
 {

   myTS = methods::new("conditionalPoissonTS", tree = as.data.frame(tree), data = as.data.frame(data), B = B)
 }
  myTS@direction = direction

  ## TBSS analysis
  myTS = mapNodesLeaves(myTS)
  myTS = computeSS(myTS)

  ## Remove nodes with a number of events lower than min_events 
   myTS@nodesToTest = myTS@nodeSS$node[(myTS@nodeSS$observed >= min_events) & (myTS@nodeSS$expected>0)] 

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

