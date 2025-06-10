#' Permutation TBSS 
#' @description This function implements a permutation based TBSS. The function uses as test the likelihood ratio test (LRT) to compare difference in mean of two Poisson distribution and can include offsets (e.g., the person time for the exposure groups). The null distribution for the LRT is computed  using permutation. Hence, unlike other TBSS this implementation require subject level data. Importantly this implementation can use node-specifyc washout, i.e. for each node it does not count patients  that have already experienced the events codified in the node in a pre-specified washout window. 
#' @param cohort a \code{matrix} or \code{data.frame} with  five or seven columns (two are optionals)  that are ('patientID', 'indexDate', 'exposure','weight', 'strata') and followupStart (0 corresponding to index data) followupEnd end of folloupTime. Strata are used for permutation and can correspond, for example, to matched observation. If all strata are equal (e.g., to 1) observations are assumed to be exchangeable. Weight is used to differentially weight observations. The followup time, if provided, is used as an offset to test differences in rates rather than means.  Note that IndexData is expected in date format (use R function as.Date() for a valid input) 
#' @param diagnosis_table a \code{matrix} or \code{data.frame} with three or four columns. The first three are 'patientID', 'date', 'leaf'. The date refer to occurrences of the event codified by leaf. Patients can have multiple events (there will be no double counting). date is expected to be in date format (use R function as.Date() for a valid input), while leaf are treated as characters. A fourth optional column ('include') is a binary variable that indicates if some events have to be censored. Censored events (time) counts in terms of person time and the events are used for the washout but not  in the analysis. 
#' @param tree a \code{matrix} or \code{data.frame} with two columns (nodeID, parentID) that specifies the hierarchical structure in an edge list format.
#' @param min_events the minimum number of observed events for a node to be tested. This is an integer greater or equal to 2. Default \code{min_events = 2}. Note that by choosing higher \code{min_events}, the power can increase because of the decreased number of hypotheses that are tested. This argument should be fixed before running the analysis to preserve the validity of the p-values.  
#' @param B number of Monte Carlo (i.e., permutations) replicates. The default is \code{B = 999}.
#' @param parallel if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.
#' @param ncpus number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.
#'@param direction can be either \code{'positive'}, \code{"negative"} or, \code{"all"} by default (\code{direction = 'positive'}) poissonTBSS only look at nodes where the number of expected cases is larger than the number of observed cases, or equivalently the excess of cases is positive. Setting  \code{direction = '"all"'} test all nodes, while \code{direction  = "negative"} only node with negative excess of cases.
#'@param nodeToRemove a character vector with user-specified nodes to be removed from the analysis.  Note that sufficient statistics for these nodes are computed, but they are excluded from the testing procedure. By removing nodes from the list of tests, power can increase because of the decreased number of hypotheses that are tested. This list should be specified before running the analysis.  
#'@param offset  indicate how expected counts are scaled to be compared with observed counts. Possible values include "sum_of_weights" indicating that expected counts are scaled with  (N weighted total exposed / N weighted total comparator) and "unscaled" for no scaling. Alternative scaling can be used by directly scaling the weights in the cohort file, and using the "unscaled" option for offset
#'@param seed a number indicating the seed used for the analysis, i.e., the input of \code{set.seed()}. Default \code{seed  = 8277}.
#' @param parsedParameters a names list of object that can be precomputed. This include mapNodeTree (i.e., the tree parsing step), and idInNode (the list of nodes/patients created by the \code{applyNodeSpecificWashout()} function).
#'@param verbose, if \code{TRUE} (default) print some messages 
#'@param relative, if \code{FALSE} force the analysis to be on means rather than rates even when folloup time start and end days are provided in the cohort file. If followup details are not provided this parameter is ignored. 

#'@export
permutationTBSS = function(cohort, diagnosis_table,tree, min_events = 2, B = 999, parallel = FALSE, ncpus = NULL, nodeToRemove = NULL,direction = 'positive', offset = "sum_of_weights", seed  = 8277, parsedParameters = NULL, verbose = TRUE, relative = TRUE)
{
match.arg(offset, c("sum_of_weights", "unscaled"))
set.seed(seed)
## stop if the parameter min_events is smaller than 2 or negative 
if(min_events < 2) stop("min_events should be a number greater than 2.")

 tree = as.data.frame(tree)

 if(parallel & (length(ncpus) == 0))
 {
  ncpus = detectCores() - 1
 }

## names data appropriatly - assumng correct order of the colum

fup = FALSE
names(cohort)[1:5] = c('patientID', 'indexDate', 'exposure','weight', 'strata')
if(NCOL(cohort)>5)
{ 
	if(NCOL(cohort)<7) stop("You need to specify both followupStart and followupEnd. The cohort data used as input do not have all the necessary input.")
	names(cohort)[6] = 'followupStart' 
	names(cohort)[7] = 'followupEnd' 
	fup = TRUE
}
cohort = cohort[order(cohort$strata),]
cohort$patientID = as.character(cohort$patientID)
cohort$indexDate = as.Date(cohort$indexDate)
## scale weights
if(offset =='sum_of_weights')
{
  wN = with(cohort,tapply(weight,exposure, sum))
  offset = wN[2]/wN[1]
  cohort$weight  = with(cohort, ifelse(exposure ==0,weight * offset, weight))    
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
diagnosis_table =  left_join(diagnosis_table, cohort, by= 'patientID') 
## create a variable indicatings date from index
diagnosis_table$time = with(diagnosis_table, date - indexDate)
if(fup)
{
 ind_input = diagnosis_table$time >= diagnosis_table$followupStart & diagnosis_table$time <= diagnosis_table$followupEnd
ind_wash = !ind_input
}
else
{
 ind_input = diagnosis_table$time >=0 
 ind_wash = !ind_input
}

inputData = diagnosis_table[ind_input & diagnosis_table$include == 1 ,]
washData = diagnosis_table[ ind_wash ,]


## appropriately scale expected events

########################
if(!relative){fup = FALSE}

myTS = new("permutationPoissonTS", tree = tree, data = data.frame(leaf = inputData$leaf, exposure = inputData$exposure ), fup = fup, B =B, cohort = cohort)

if(verbose) cat("Parsing the tree...")
tim = system.time({
if('mapNodesLeaves' %in% names(parsedParameters))
{
	myTS@mapNodesLeaves = parsedParameters$mapNodesLeaves
}
else{myTS = TBSS:::mapNodesLeaves(myTS)}
})
if(verbose) cat(sprintf("Done. \n \n Computation took approximately %.3f seconds \n\n\n", tim[3]))

if(verbose) cat("Applying node-specific washout, this can take some time...")

tim = system.time({
if('idInNode' %in% names(parsedParameters))
{
myTS@idInNode = parsedParameters$idInNode
}
else
{
 myTS = applyNodeSpecificWashout(myTS,washData, inputData)
}
})
if(verbose) cat(sprintf("Done. \n\n Computation took approximately %.3f seconds \n\n\n", tim[3]))

myTS = computeSS(myTS) 
## Remove nodes with a number of events lower than min_events 
   if(direction == 'positive')
   {
    myTS@nodesToTest = myTS@nodeSS$node[(myTS@nodeSS$observed >= min_events) & ((myTS@nodeSS$observed - myTS@nodeSS$expected)>0) & (myTS@nodeSS$expected>0)]
   }
   else if (direction == 'negative')
   {
     myTS@nodesToTest = myTS@nodeSS$node[(myTS@nodeSS$observed >= min_events) & ((myTS@nodeSS$observed - myTS@nodeSS$expected)<0 &(myTS@nodeSS$expected>0) )] 
   }
   else
   {
     myTS@nodesToTest = myTS@nodeSS$node[(myTS@nodeSS$observed >= min_events) & (myTS@nodeSS$expected>0)] 
   }

  ## Remove user specified node
  if(!is.null(nodeToRemove))
  {
   myTS@nodesToTest = setdiff(myTS@nodesToTest,nodeToRemove) 
  }

myTS = computeLRT(myTS) 



## perform MC

if(verbose) cat("Starting permutation algorithm...")
tim = system.time({myTS = monteCarlo(myTS, parallel, ncpus)})
if(verbose) cat(sprintf("Done. \n\n Computation took approximately %.3f seconds \n\n\n", tim[3]))
## return results
return(myTS)
}

