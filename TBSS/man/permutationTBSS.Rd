% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/permutationTBSS.R
\name{permutationTBSS}
\alias{permutationTBSS}
\title{Permutation TBSS}
\usage{
permutationTBSS(
  cohort,
  diagnosis_table,
  tree,
  min_events = 2,
  B = 999,
  parallel = FALSE,
  ncpus = NULL,
  nodeToRemove = NULL,
  direction = "positive",
  offset = "sum_of_weights",
  seed = 8277,
  parsedParameters = NULL,
  verbose = TRUE,
  relative = TRUE
)
}
\arguments{
\item{cohort}{a \code{matrix} or \code{data.frame} with  five or seven columns (two are optionals)  that are ('patientID', 'indexDate', 'exposure','weight', 'strata') and followupStart (0 corresponding to index data) followupEnd end of folloupTime. Strata are used for permutation and can correspond, for example, to matched observation. If all strata are equal (e.g., to 1) observations are assumed to be exchangeable. Weight is used to differentially weight observations. The followup time, if provided, is used as an offset to test differences in rates rather than means.  Note that IndexData is expected in date format (use R function as.Date() for a valid input)}

\item{diagnosis_table}{a \code{matrix} or \code{data.frame} with three or four columns. The first three are 'patientID', 'date', 'leaf'. The date refer to occurrences of the event codified by leaf. Patients can have multiple events (there will be no double counting). date is expected to be in date format (use R function as.Date() for a valid input), while leaf are treated as characters. A fourth optional column ('include') is a binary variable that indicates if some events have to be censored. Censored events (time) counts in terms of person time and the events are used for the washout but not  in the analysis.}

\item{tree}{a \code{matrix} or \code{data.frame} with two columns (nodeID, parentID) that specifies the hierarchical structure in an edge list format.}

\item{min_events}{the minimum number of observed events for a node to be tested. This is an integer greater or equal to 2. Default \code{min_events = 2}. Note that by choosing higher \code{min_events}, the power can increase because of the decreased number of hypotheses that are tested. This argument should be fixed before running the analysis to preserve the validity of the p-values.}

\item{B}{number of Monte Carlo (i.e., permutations) replicates. The default is \code{B = 999}.}

\item{parallel}{if \code{parallel = TRUE} parallel computing is used to compute the Monte Carlo p-values.}

\item{ncpus}{number of cpus used if \code{parallel = TRUE}. If no value is provided, the default is \code{ncpus = detectCores() - 1}.}

\item{nodeToRemove}{a character vector with user-specified nodes to be removed from the analysis.  Note that sufficient statistics for these nodes are computed, but they are excluded from the testing procedure. By removing nodes from the list of tests, power can increase because of the decreased number of hypotheses that are tested. This list should be specified before running the analysis.}

\item{direction}{can be either \code{'positive'}, \code{"negative"} or, \code{"all"} by default (\code{direction = 'positive'}) poissonTBSS only look at nodes where the number of expected cases is larger than the number of observed cases, or equivalently the excess of cases is positive. Setting  \code{direction = '"all"'} test all nodes, while \code{direction  = "negative"} only node with negative excess of cases.}

\item{offset}{indicate how expected counts are scaled to be compared with observed counts. Possible values include "sum_of_weights" indicating that expected counts are scaled with  (N weighted total exposed / N weighted total comparator) and "unscaled" for no scaling. Alternative scaling can be used by directly scaling the weights in the cohort file, and using the "unscaled" option for offset}

\item{seed}{a number indicating the seed used for the analysis, i.e., the input of \code{set.seed()}. Default \code{seed  = 8277}.}

\item{parsedParameters}{a names list of object that can be precomputed. This include mapNodeTree (i.e., the tree parsing step), and idInNode (the list of nodes/patients created by the \code{applyNodeSpecificWashout()} function).}

\item{verbose, }{if \code{TRUE} (default) print some messages}

\item{relative, }{if \code{FALSE} force the analysis to be on means rather than rates even when folloup time start and end days are provided in the cohort file. If followup details are not provided this parameter is ignored.}
}
\description{
This function implements a permutation based TBSS. The function uses as test the likelihood ratio test (LRT) to compare difference in mean of two Poisson distribution and can include offsets (e.g., the person time for the exposure groups). The null distribution for the LRT is computed  using permutation. Hence, unlike other TBSS this implementation require subject level data. Importantly this implementation can use node-specifyc washout, i.e. for each node it does not count patients  that have already experienced the events codified in the node in a pre-specified washout window.
}
