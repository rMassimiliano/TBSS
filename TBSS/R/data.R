#' An example of tree with 25 nodes
#'
#' Describes a hierarchical tree structure used as standard input in TBSS.
#'
#' @format A data frame with 24 rows and 2 columns describing a hierarchical structure:
#' \describe{
#'   \item{childID}{node ID from where an edge starts, character.}
#'   \item{parentID}{node ID to where an edge end, character.}
#' }
#' @source Simulated dataset.
"tree_example"

#' Example dataset for Bernoulli TBSS
#' A list with two elements. Each element contains the case and control file for implementation of a Bernoulli TBSS. The first element of the list corresponds to the standard Bernoulli TBSS while the second consider variable matching ratio. 
#' @format A named list containing two data frames:
#' \describe{
#'   \item{fixed_matching_ratio}{
#'     Two data frame with 13 rows, one for the case on for the control and the probability (p) used in the example.
#'     \describe{
#'       \item{leaf}{Unique identifier for the leaf, character}
#'       \item{case}{number of cases, integer}
#'       \item{control}{number of controls, integer}
#'     }
#'   }
#'   \item{variable_matching_ratio}{
#'	Two data frames with 37 rows and 3 columns 
#'     \describe{
#'       \item{leaf}{Unique identifier for the leaf, character}
#'       \item{case}{number of cases, integer}
#'       \item{control}{number of controls, integer}
#'       \item{p}{probability associated with corresponding leaf and strata, numeric}
#'     }
#'   }
#' }
#' @source Simulated dataset.
"bernoulli_TBSS_examples"


#' Example dataset for Poisson TBSS
#' @format A data frame with 13 rows and 3 columns
#'     \describe{
#'       \item{leaf}{Unique identifier for the leaf, character}
#'       \item{observed}{number of observed cases, integer}
#'       \item{expected}{number of expected cases (appropriately scaled), numeric}
#'     }
#' @source Simulated dataset.
"poisson_TBSS_example"

#' Example dataset for Poisson (exponential) and robust survival TBSS
#' @format The format of the data is a data frame with 13 rows and 11 columns
#'     \describe{
#'       \item{leaf}{Unique identifier for the leaf, character}
#'	 \item{sample_size_case}{number of individuals for the cases, numeric}
#'	 \item{events_case}{number of events for cases, numeric}
#'	 \item{total_person_time_case}{total person time for the cases, numeric}
#'	 \item{total_person_time_SOS_case}{sum of squares of the total person time for the cases, numeric}
#'	 \item{uncensored_person_time_case}{total person time among the uncensored individuals for the cases, numeric}
#'	 \item{sample_size_control}{number of individual for the controls, numeric}
#'	 \item{events_control}{number of events for control, numeric}
#'	 \item{total_person_time_control}{total person time for the controls, numeric}
#'	 \item{total_person_time_SOS_control}{sum of squares of the total person time for the controls, numeric}
#'	 \item{uncensored_person_time_control}{total person time among the uncensored individuals for the controls, numeric}
#'     }
#' @source Simulated dataset.
"survival_poisson_example"

#' Example dataset for permutation based TBSS (permutationPoissonTS) 
#' @format The format of the data is a list with two data sets cohort and diagnosis table
#'     \describe{
#'       \item{cohort}{ a \code{data.frame} with seven columns 
#'	 that are 'patientID', 'indexDate', 'exposure','weight', 'strata', 'followupStart' and 'followupEnd'.}  
#'	 \item{diagnosis table}{ a \code{data.frame} with three columns: 'patientID', 'date', 'leaf'. The date refer to occurrences of the event codified by leaf. Patients can have multiple events (there will be no double counting). }
#'}
#' @source Simulated dataset.
"patient_level_data"



#' Example for Case-Case-Time-Control

#' @format A named list containing two data frames:
#' \describe{
#'     Two data frame with 13 rows, one for the case on for the control 
#'     \describe{
#'       \item{leaf}{Unique identifier for the leaf, character}
#'	\item{control_reference}{Number of events for the controls in the reference window}
#' 	\item{control_hazard}{Number of events for the controls in the hazard window}
#'	\item{case_reference}{Number of events for the cases in the reference window}
#' 	\item{case_hazard}{Number of events for the cases in the hazard window}
#'     }
#'   }
#' @source Simulated dataset.
"cctc_TBSS_example"
