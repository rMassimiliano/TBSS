#' Node specific washout
#' @description function that applies a node specific washout. 
#' @param object a \code{permutationPoissonTS} or other appropriate object derived from the class TS. 
#' @param washData a \code{data.frame} with the unique combination of patient IDs and codes at the leaf level of the tree
#'@param inputData data including events happened from the index data until end of follow-up
#'@details This function assumes that some slots of the TS class have been appropriately initialized. Specifically the: 1) \code{cohort} slot where each row is a patient, columns are (patientID,  indexDate, exposure,   weight, follow-up); 2) the \code{mapNodeLeaves} slot that is a list that for each node provides the corresponding leaves linked to that node; 3) the  \code{fup} slot that if TRUE compute the person time for each individual in each node.
#'@returns The function returns as output an updated TS object.
#' The updated TS object includes a the slot \code{idInNode} which is a named list where each element is a tree node associated with a vector o patientsID to be considered in that node
#' If \code{object@fup = TRUE}, that is to say we are accounting for person time, other two named list \code{personTime} and \code{personTimeNotInNode} are also included in the output. 
#' The \code{personTime} includes the total person time, in days, for the exposed and comparator groups.
#'The person time is weighted using the weight provided in the cohort file. In case an event happens the same day as the index date, we assume that the event after 1/2 a day. 
#' Patients excluded from the washout procedure are "censored" rather than excluded from the analysis, hence their person-time counts. For "censored" patients, we censor their time at the first occurring event or the length of follow-up if they experience no events.
#'\code{personTimeNotInNode} includes person time for all the individuals that did not experienced that node event.
#' @references
#'Russo, M., Sreedhara, S.K., Smith, J., Davis, S.E., Maro, J.C., Deramus, T., Lii, J., Yang, J., Desai, R.J., Hernández-Muñoz, J.J. and Ma, Y., 2026. Electronic health record-enhanced signal detection using tree-based scan statistic methods. American Journal of Epidemiology, 195(1), pp.178-187.
#' @importFrom rlang .data
applyNodeSpecificWashout = function(object,washData, inputData)
{
 retval = list()
 personTime = list()
 personTimeNotInR = list()
for(r in 1:length(object@mapNodesLeaves))
{
  ## unique patients IDs to keep 
   idToRm = unique(washData[with(washData, leaf %in% object@mapNodesLeaves[[r]]),]$patientID)
  ## All patient  IDs 
   allNodeID = unique(inputData[with(inputData, leaf %in% object@mapNodesLeaves[[r]]),]$patientID)
  if(NROW(idToRm)>0)
  {
   ## all events are the same
   retval[[r]] = object@cohort[object@cohort$patientID %in% setdiff(allNodeID,idToRm),]
  }
  else
  {
   retval[[r]] = object@cohort[object@cohort$patientID %in% allNodeID,]
  }
  if(object@fup)
  {
    if(NROW(retval[[r]])>0)
    {
        ## all events for a certain node from the followup time to end
        tmp =  inputData[inputData$patientID %in% retval[[r]]$patientID,]
        tmp = tmp[tmp$leaf %in% object@mapNodesLeaves[[r]],]
        ## use shorter event time certain node
        fups = tapply(tmp$time - tmp$followupStart, tmp$patientID,min)
        fups = data.frame(patientID = names(fups), followup = fups)
        ## happening the same day event time 0.5 days
        fups$followup[fups$followup==0] = 0.5
    }
    else
    {
        fups = data.frame(patientID = character(0), followup = numeric(0))
    }
  
    ## all events for "censored" patients
    tmp2 = inputData[inputData$patientID %in% setdiff(allNodeID,retval[[r]]$patientID),]
    tmp2 = tmp2[tmp2$leaf %in% object@mapNodesLeaves[[r]],]

    if(NROW(tmp2)>0)
    {
      ## followup-for censored individuals
      fups2 = tapply(tmp2$time - tmp2$followupStart,tmp2$patientID,min)
      fups2 = data.frame(patientID = names(fups2), followup = fups2)
    }
    else
    {
	    ## copy the structure in a null data.frame
	    fups2 = fups[0,]
    }

    ## person time (full follow-up length) for individuals not associated with node r at all
    notInR = object@cohort[!(object@cohort$patientID %in% allNodeID),]
    notInR$followup = with(notInR, followupEnd - followupStart)
    notInR$wft = with(notInR, weight * followup)
    personTimeNotInR[[r]] = notInR[, c('patientID','followup','exposure','weight','wft')]

    # add folloup time to the "sufficient" statistics -- no need  because the information is in person time
    retval[[r]]  =  dplyr::left_join(retval[[r]], fups, by = 'patientID')



   ptime = dplyr::left_join(fups2, object@cohort |> dplyr::filter(.data$patientID %in% fups2$patientID) |> dplyr::select(.data$patientID, .data$exposure,.data$weight),  by = 'patientID')

   ptime$wft = with(ptime,weight * followup)

   #personTime[[r]] = tapply(ptime$wft,ptime$exposure,sum)
   personTime[[r]] = ptime
  }

}
object@idInNode = retval
names(object@idInNode) = names(object@mapNodesLeaves)
if(object@fup)
{
 object@personTime = personTime
 names(object@personTime) = names(object@mapNodesLeaves)
 object@personTimeNotInNode = personTimeNotInR
 names(object@personTimeNotInNode) = names(object@mapNodesLeaves)
}

return(object)
}
