## This file pre-define the classes used in the package to avoid issues due to the fact that files are loaded in alphabetical order
#General TS class for tree scan
#'General TS class for TBSS
#' General TS class.
#' @description This class is used as a template for all the TBSS implementation. Several generic function (e.g. utilities to parse the tree) are implemented in this generic class
#'@slot data contains the data in a \code{data.frame} one row for each leaf of the tree.  The number of columns varies depending on the analysis, while the first column is the name of a leaf. 
#'@slot tree a\code{data.frame} with two columns representing the hierarchical structure in a edge list format. Each row represents a directional edge between a node and its parent node, so the two columns are a nodeID and parentID. The IDs are treated as \code{character}.
#'@slot leaves is a vector with the names of the leaves, typically the first column of the slot \code{data}
#'@slot nodes a vector with the names of all the nodes in the tree, i.e. the first column of \code{tree}
#'@slot mapNodesLeaves is a named list. Each element of the list is named after a node, while the content of the list indicates the leaves associated with that node. 
#'@slot nodesToTest is a vector defining the nodes to be tested. Some nodes could be excluded from the testing procedure  (e.g., because there are too few data points) and can be excluded from this list.
#'@slot nodeSS is a \code{data.frame} reporting the sufficient statistics for each node in the tree. The number of columns varies depending on the analysis, and the first column is the name of a node. 
#'@slot LRT is a named vector containing the node-specific discrepancy measure computed on the observed data for each node in \code{nodesToTest}.
#' @slot B is the number of MC replicates used to approximate the null distribution of LRT
#'@slot LRT_H0 is a \code{matrix} with the MC-distribution of LRT. Each column is an MC replicate.
#' @docType class
#' @name TS-class
#' @aliases TS-class TS
setClass(Class = "TS", slots = c(data = 'data.frame',
                                 tree ='data.frame',
                                 leaves = 'character',
                                 nodes = 'character',
                                 mapNodesLeaves = 'list',
				 nodesToTest ='character',
                                 nodeSS = 'data.frame',
                                 LRT ='numeric',
                                 B = 'numeric',
                                 LRT_H0 = 'matrix'))

#' Unconditional Bernoulli TS class
#' @description Define the class for Unconditional Bernoulli TBSS and directly inherit from TS
#' @slot p represent the probability of being case/control used for the analysis
#' @docType class
#' @name unconditionalBernoulliTS
setClass(Class = "unconditionalBernoulliTS", slots = c(p = 'numeric'), contains = 'TS')

#' Unconditional Bernoulli with Variable matching ratio TS class
#' @description Define the class for Unconditional Bernoulli TBSS and directly inherit from TS
#' @docType class
#' @name unconditionalBernoulliTSVariableRR
setClass(Class = "unconditionalBernoulliTSVariableRR", contains = 'unconditionalBernoulliTS')

#' Unconditional Case Case time Control
#' @description Define the class for Unconditional Case Case time Control Bernoulli TBSS and directly inherit from TS
#' @slot p represent the probability of being case/control used for the analysis
#' @docType class
#' @name unconditionalCaseCaseTimeControlTS
setClass(Class = "unconditionalCaseCaseTimeControlTS", slots = c(p = 'numeric'), contains = 'TS')


#' Unconditional Poisson TS class
#' @description Define the class for Unconditional Poisson TBSS and directly inherit from TS
#' @docType class
#' @name unconditionalPoissonTS
setClass(Class = "unconditionalPoissonTS",  contains = 'TS')

#' Conditional Poisson TS class
#' @description Define the class for Unconditional Poisson TBSS and inherit from unconditionalPoissonTS
#' @docType class
#' @name conditionalPoissonTS 
setClass(Class = "conditionalPoissonTS",  contains = 'unconditionalPoissonTS')


#'  Permutation TBSS 
#' @description This class is used to implement permutation TBSS with node specific washout 
#'@slot data contains the data in a \code{data.frame} one row for each leaf of the tree.  The number of columns varies depending on the analysis, while the first column is the name of a leaf. 
#'@slot tree a\code{data.frame} with two columns representing the hierarchical structure in a edge list format. Each row represents a directional edge between a node and its parent node, so the two columns are a nodeID and parentID. The IDs are treated as \code{character}.
#'@slot idInNode  list information of individuals in each node of the tree. This is internally created by the function `applyNodeSpecificWashout()`
#'@slot personTime list personTime of individuals in each node of the tree whose events are washed out. This is internally created by the function `applyNodeSpecificWashout()`
#'@slot personTimeNotInNode list personTime of individuals NOT is an node. This is internally created by the function `applyNodeSpecificWashout()`
#' @param cohort a \code{matrix} or \code{data.frame} with  five or seven columns (two are optionals)  that are ('patientID', 'indexDate', 'exposure','weight', 'strata') and followupStart (0 corresponding to index data) followupEnd end of folloupTime. Strata are used for permutation and can correspond, for example, to matched observation. If all strata are equal (e.g., to 1) observations are assumed to be exchangeable. Weight is used to differentially weight observations. The followup time, if provided, is used as an offset to test differences in rates rather than means.  Note that IndexDate is expected in date format (use R function as.Date() for a valid input). 
#'@slot fup indicate if time is used in the analysis. If false only number of events are used.
#' @docType class
#' @name permutationPoissonTS-class
#' @aliases permutationPoissonTS-class permutationPoissonTS
#' @importFrom rlang .data
setClass(Class = "permutationPoissonTS",  contains = 'TS', slots = c(idInNode = 'list',personTime = 'list',personTimeNotInNode = 'list', cohort = 'data.frame', fup = 'logical')) 
