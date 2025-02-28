#' This function implement TreeScan for a case-case-time control design
#' @param `case`: a matrix or data.frame with three columns (leavesID, cases in the reference window, cases in the hazard window)
#' @param `control`: a matrix or data.frame with three columns (leavesID, cases in the reference window, cases in the hazard window)
#' @param `tree`: a matrix or data.frame with two columns (nodeID, parantNodeID) 
#' @param `B` number of bootstrap replicates. Default is `999`.
#' @param `parallel` if TRUE parallel computing is used to compute the bootstrap p-values.
#' @param `ncpus` number of cpus used if `parallel = TRUE`. If no value is provided set to  `detectCores() - 1`.

#' @export
cctcTBSS =  function(case,control,tree,B = 999, parallel = FALSE, ncpus = NULL, nodeToRemove = NULL, seed  = 8277)
{
 ## create data 
 names(case) = c("leaf", "case_reference", "case_hazard")
 names(control) = c("leaf", "control_reference", "control_hazard")
 data = control |> dplyr::left_join(case, by = 'leaf') |> as.data.frame()

 tree = as.data.frame(tree)
 if(parallel & (length(ncpus) == 0))
 {
  ncpus = detectCores() - 1
 }

 ## create appropriate TS object
  myTS = new("unconditionalCaseCaseTimeControlTS", data =data, tree = tree, B =B )
  myTS = mapNodesLeaves(myTS)
  myTS = computeSS(myTS)
  myTS = computeLRT(myTS)
  myTS = monteCarlo(myTS, parallel, ncpus)

  return(myTS)
}

