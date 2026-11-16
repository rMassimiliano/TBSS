#' Permutation TBSS 
#' @description This function implements a permutation based TBSS. The function uses as test the likelihood ratio test (LRT) to compare difference in mean of two Poisson distribution and can include offsets (e.g., the person time for the exposure groups). The null distribution for the LRT is computed  using permutation. Hence, unlike other TBSS this implementation require subject level data. Importantly this implementation can use node-specifyc washout, i.e. for each node it does not count patients  that have already experienced the events codified in the node in a pre-specified washout window. 
#' @param cohort a \code{matrix} or \code{data.frame} with  five or seven columns (two are optionals)  that are ('patientID', 'indexDate', 'exposure','weight', 'strata') and followupStart (0 corresponding to index data) followupEnd end of folloupTime. Strata are used for permutation and can correspond, for example, to matched observation. If all strata are equal (e.g., to 1) observations are assumed to be exchangeable. Weight is used to differentially weight observations. The followup time, if provided, is used as an offset to test differences in rates rather than means.  Note that IndexDate is expected in date format (use R function as.Date() for a valid input). 
#' @param diagnosis_table a \code{matrix} or \code{data.frame} with three or four columns. The first three are 'patientID', 'date', 'leaf'. The date refer to occurrences of the event codified by leaf. Patients can have multiple events (there will be no double counting). date is expected to be in date format (use R function as.Date() for a valid input), while leaf are treated as characters. A fourth optional column ('include') is a binary variable that indicates if some events have to be censored. Censored events (time) counts in terms of person time and the events are used for the washout but not  in the analysis. 
#' @param tree a \code{matrix} or \code{data.frame} with two columns (nodeID, parentID) that specifies the hierarchical structure in an edge list format.
#' @param min_events the minimum number of observed events for a node to be tested. This is an integer greater or equal to 2. Default \code{min_events = 2}. Note that by choosing higher \code{min_events}, the power can increase because of the decreased number of hypotheses that are tested. This argument should be fixed before running the analysis to preserve the validity of the p-values.  
#' @param B number of Monte Carlo (i.e., permutations) replicates. The default is \code{B = 999}.
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.
#'@param direction can be either \code{'positive'}, \code{"negative"} or \code{"all"}. By default (\code{direction = 'positive'}), the LRT of the nodes where the observed number of cases (or rate) is not larger than the expected one is set to 0. The restriction is applied in the same way to the observed data and to every Monte Carlo replicate. Setting \code{direction = "all"} tests all nodes, while \code{direction = "negative"} only considers nodes where the observed number of cases (or rate) is not smaller than the expected
#'@param nodeToRemove a character vector with user-specified nodes to be removed from the analysis.  Note that sufficient statistics for these nodes are computed, but they are excluded from the testing procedure. By removing nodes from the list of tests, power can increase because of the decreased number of hypotheses that are tested. This list should be specified before running the analysis.  
#'@param washout_window  a numeric vector with two negative numbers that indicate the length  of the washout window relative to the index data. For example, if washout_window = c(-180,-1) all diagnosis events that occurred in the 180 days (including days -180 and -1)  before the index are used for washout. The default is \code{washout_window = NULL} indicating that all events before the index date (excluded) are used for washout.
#'@param offset  indicate how expected counts are scaled to be compared with observed counts. Possible values include "sum_of_weights" indicating that expected counts are scaled with  (N weighted total exposed / N weighted total comparator) and "unscaled" for no scaling. Alternative scaling can be used by directly scaling the weights in the cohort file, and using the "unscaled" option for offset. By default no scaling is applied.
#' @param statistic indicate the test statistic to be used for the permutation test. Possible values are `Poisson` (the default) that uses the log-likelihood ratio test for testing the difference of the ratios of two Poisson distributions and `Cox` that uses the log partial likelihood ratio from a Cox model. 
#'@param seed a number indicating the seed used for the analysis, i.e., the input of \code{set.seed()}. Default \code{seed  = 8277}.
#' @param parsedParameters a names list of object that can be precomputed. This include mapNodeTree (i.e., the tree parsing step), and idInNode (the list of nodes/patients created by the \code{applyNodeSpecificWashout()} function). When follow-up time is specified, personTime and personTimeNotInNode are also needed. These parameters are created by the function \code{applyNodeSpecificWashout()} and if they are not provided they are recomputed.
#'@param verbose, if \code{TRUE} (default) print some messages 
#'@param relative, if \code{FALSE} force the analysis to be on means rather than rates even when followup time start and end days are provided in the cohort file. If followup details are not provided this parameter is ignored. 
#'@details Note if data regarding followup time is not provided we use the convention that all events in the diagnosis table occurring from the index date on (with the index date excluded) are used for the analysis. Similarly, all data preceding the index date (included) are used for washout. This behavior can be changed by setting the \code{washout_window} argument and setting  the columns corresponding to the parameters for the follow-up in the cohort file appropriately.
#'@references 
#' Russo, M., Sreedhara, S.K., Smith, J., Davis, S.E., Maro, J.C., Deramus, T., Lii, J., Yang, J., Desai, R.J., Hernández-Muñoz, J.J. and Ma, Y., 2026. Electronic health record-enhanced signal detection using tree-based scan statistic methods. American Journal of Epidemiology, 195(1), pp.178-187.
#'
#' Russo, M., Hahn, G., Huybrechts, K. F., & Wang, S. V. (2026). A tree-based scan statistic for database studies with time-to-event outcomes. arXiv preprint arXiv:2609.38582.
#' @export
#' @examples
#' library(TBSS)
## load data
#' data(patient_level_data, package = 'TBSS')
#' attach(patient_level_data)
#' ## load tree
#' data(tree_example, package = 'TBSS')
#' 
#' ## running the algorithm eith 30 permutation with 2 cores in parallel for this example
#' ## default is 999 without parallelization
#' mod_tbss = permutationTBSS(cohort, diagnosis_table, tree_example,
#'			       B = 30, parallel = TRUE, ncpus = 2) 
#' 
#' summary(mod_tbss) |> head()
permutationTBSS = function(cohort, diagnosis_table,tree, min_events = 2, B = 999, parallel = FALSE, ncpus = NULL, nodeToRemove = NULL,washout_window = NULL, direction = 'positive', offset = "unscaled", statistic = 'Poisson', seed  = 8277, parsedParameters = NULL, verbose = TRUE, relative = TRUE)
{
set.seed(seed)
offset = match.arg(offset, c("sum_of_weights", "unscaled"))
statistic = match.arg(statistic, c("Poisson", "Cox"))
## stop if the direction is not among the possibilities
if(!(direction %in% c("positive","all","negative"))) stop("direction should be 'positive', 'all', or 'negative'. Refer to the help for details.")

## stop if the parameter min_events is smaller than 2 or negative 
if(min_events < 2) stop("min_events should be a number greater than or equal to 2.")

 tree = as.data.frame(tree)

 if(parallel & (length(ncpus) == 0))
 {
  ncpus = parallel::detectCores() - 1
 }

## names data appropriately - assuming correct order of the columns

fup = FALSE
names(cohort)[1:5] = c('patientID', 'indexDate', 'exposure','weight', 'strata')
if(NCOL(cohort)>5)
{ 
	if(NCOL(cohort)<7) stop("You need to specify both followupStart and followupEnd. The cohort data used as input do not have all the necessary input.")
	names(cohort)[6] = 'followupStart' 
	names(cohort)[7] = 'followupEnd' 
	fup = TRUE

## check that there are no inconsistencies in the definition of followup
if(anyNA(cohort$followupStart) | anyNA(cohort$followupEnd)) stop("followupStart and followupEnd cannot contain missing values.")
if(any(cohort$followupEnd < cohort$followupStart)) stop("followupEnd should not be smaller than followupStart  (no negative length of the follow-up).")
if(!is.null(washout_window) && any(washout_window[2] >= cohort$followupStart)) warning("washout_window overlaps the follow-up of some patients: their events in the overlap are used both for the washout and in the analysis.")
}


if(statistic == 'Cox' & !fup) stop("statistic = 'Cox' requires follow-up time: provide followupStart and followupEnd in the cohort data.")
cohort$patientID = as.character(cohort$patientID)
cohort$indexDate = as.Date(cohort$indexDate)
## scale weights
if(offset =='sum_of_weights')
{
  wN = with(cohort,tapply(weight,exposure, sum))
  offset = wN[2]/wN[1]
  cohort$weight  = with(cohort, ifelse(exposure ==0,weight * offset, weight))    
  if(statistic == 'Cox') warning("Are you sure you want to scale weights while fitting a Cox model? consdier using offset == 'unscaled'")
}

names(diagnosis_table)[1:3] = c('patientID', 'date', 'leaf')
if(NCOL(diagnosis_table) > 3)
{
 names(diagnosis_table)[4] = "include"
}
else
{
    diagnosis_table$include = 1
}


## "convert" data in appropriate format
diagnosis_table$date = as.Date(diagnosis_table$date)
diagnosis_table$patientID = as.character(diagnosis_table$patientID)
diagnosis_table$leaf = as.character(diagnosis_table$leaf)


## remove from the diagnosis table patients that are not included in the cohort (if any)
to_keep = which((diagnosis_table$patientID %in% cohort$patientID))
diagnosis_table = diagnosis_table[to_keep,]

## 2nd divide diagnosis table in two: pre and post washout window
 # add the index date from cohort
diagnosis_table =   dplyr::left_join(diagnosis_table, cohort, by= 'patientID') 
## create a variable indicating date from index
diagnosis_table$time = with(diagnosis_table, date - indexDate)
if(fup)
{
 ind_input = diagnosis_table$time >= diagnosis_table$followupStart & diagnosis_table$time <= diagnosis_table$followupEnd
 ## if washout_window is NULL then we use all events before the index date
 ind_wash = if(is.null(washout_window)) diagnosis_table$time < diagnosis_table$followupStart  else diagnosis_table$time >= washout_window[1] & diagnosis_table$time <= washout_window[2]
 #ind_wash = !ind_input
}
else
{
 ind_input = diagnosis_table$time >0 
 ind_wash = if(is.null(washout_window)) diagnosis_table$time <=0 else diagnosis_table$time >= washout_window[1] & diagnosis_table$time <= washout_window[2]
}

inputData = diagnosis_table[ind_input & diagnosis_table$include == 1 ,]
washData = diagnosis_table[ind_wash ,]


## appropriately scale expected events

########################
if((!relative) & (statistic != 'Cox')){fup = FALSE}

if(statistic == 'Poisson')
{
myTS = methods::new("permutationPoissonTS", tree = tree, data = data.frame(leaf = inputData$leaf, exposure = inputData$exposure ), fup = fup, B =B, cohort = cohort)
}
else if(statistic == 'Cox')
{
myTS = methods::new("coxSurvivalTS", tree = tree, data = data.frame(leaf = inputData$leaf, exposure = inputData$exposure ), fup = fup, B =B, cohort = cohort)
}
myTS@direction = direction
if(verbose) cat("Parsing the tree...")
tim = system.time({
if('mapNodesLeaves' %in% names(parsedParameters))
{
	myTS@mapNodesLeaves = parsedParameters$mapNodesLeaves
}
else{myTS = mapNodesLeaves(myTS)}
})
if(verbose) cat(sprintf("Done. \n \n Computation took approximately %.3f seconds \n\n\n", tim[3]))

if(verbose) cat("Applying node-specific washout, this can take some time...")

tim = system.time({
needed = if(fup) c('idInNode', 'personTime', 'personTimeNotInNode') else 'idInNode'
if(all(needed %in% names(parsedParameters)))
{
   myTS@idInNode = parsedParameters$idInNode
   if(fup)
   {
    myTS@personTime = parsedParameters$personTime
    myTS@personTimeNotInNode = parsedParameters$personTimeNotInNode
   }
   lists = if(fup) list(myTS@idInNode, myTS@personTime, myTS@personTimeNotInNode) else list(myTS@idInNode)
   if(!all(sapply(lists, \(x) all(names(myTS@mapNodesLeaves) %in% names(x))))) stop("parsedParameters: idInNode, personTime and personTimeNotInNode should be lists named by the nodes of the tree, as returned by applyNodeSpecificWashout().")
}
else
{
 myTS = applyNodeSpecificWashout(myTS,washData, inputData)
}
})
if(verbose) cat(sprintf("Done. \n\n Computation took approximately %.3f seconds \n\n\n", tim[3]))

myTS = computeSS(myTS) 
## Remove nodes with a number of events lower than min_events 
myTS@nodesToTest = myTS@nodeSS$node[(myTS@nodeSS$observed >= min_events) & (myTS@nodeSS$expected>0)] 

## Remove user specified node
if(!is.null(nodeToRemove))
{
 myTS@nodesToTest = setdiff(myTS@nodesToTest,nodeToRemove) 
}
if(length(myTS@nodesToTest) == 0) stop("No outcome node to be tested after removing nodes from the list for min_events conditions and explicitely with nodeToRemove.")

myTS = computeLRT(myTS) 

## perform MC
if(verbose) cat("Starting permutation algorithm...")
tim = system.time({myTS = monteCarlo(myTS, parallel, ncpus)})
if(verbose) cat(sprintf("Done. \n\n Computation took approximately %.3f seconds \n\n\n", tim[3]))
## return results
return(myTS)
}

