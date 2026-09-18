
#' An example of tree with 25 nodes
#'
#' The format of the data is the standard format used by TBSS with two columns, childID and parentID.   
#' 
#'
#' @format ## `tree_example`
#' A data frame with 24 rows and 2 columns describing a hierarchical structure
#' \describe{
#'   \item{childID}{node ID from where an edge start, character.}
#'   \item{parentID}{node ID to where an edge end, character.}
#' }
"tree_example"

#' Example dataset for Bernoulli TBSS
#' The file is a list with two elements. Each elements contains the case and control file for implementation of a Bernoulli TBSS. The first element of the list correspond to the standard Bernoulli TBSS while the second consider variable matching ratio. 
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
#'	Two dataframes with 37 rows and 3 columns 
#'     \describe{
#'       \item{leaf}{Unique identifier for the leaf, character}
#'       \item{case}{number of cases, integer}
#'       \item{control}{number of controls, integer}
#'       \item{p}{probability associated with corresponding leaf and strata, numeric}
#'     }
#'   }
#' }
"bernoulli_TBSS_examples"


#' Example dataset for Poisson TBSS
#' The format of the data is a dataframe with 13 rows and 3 columns
#'     \describe{
#'       \item{leaf}{Unique identifier for the leaf, character}
#'       \item{observed}{number of observed cases, integer}
#'       \item{control}{number of expected cases (appropriately scaled), numeric}
#'     }
"poisson_TBSS_example"

#' Example dataset for Poisson (exponential) and robust survival TBSS
#' The format of the data is a dataframe with 13 rows and 11 columns
#'     \describe{
#'       \item{leaf}{Unique identifier for the leaf, character}
#'	 \item{sample_size_case}{number of individual for the cases, numeric}
#'	 \item{events_case}{number of events for cases, numeric}
#'	 \item{total_person_time_case}{total person time for the cases, numeric}
#'	 \item{total_person_time_SOS_case}{sum of squares of the total person time for the cases, numeric}
#'	 \item{uncensored_person_time_case}{total person time among the uncensored individuals for the cases, numeric}
#'	 \item{sample_size_control }{number of individual for the controls, numeric}
#'	 \item{events_control}{number of events for cases, numeric}
#'	 \item{total_person_time_control}{total person time for the controls, numeric}
#'	 \item{total_person_time_SOS_control}{sum of squares of the total person time for the controls, numeric}
#'	 \item{uncensored_person_time_control}{total person time among the uncensored individuals for the controls, numeric}
#'     }
"survival_poisson_example"
