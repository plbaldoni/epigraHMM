#' Perform peak calling of epigenomic data sets
#'
#' This function runs either consensus (one condition, multiple samples) or
#' differential (multiple conditions and samples) peak callers for epigenomic
#' data.
#'
#' @param object an epigraHMMDataSet
#' @param control list of control arguments from \code{\link{controlEM}}
#' @param type character, either \code{"consensus"} or \code{"differential"}
#' @param dist character, either \code{"zinb"} or \code{"nb"} (default)
#'
#' @return
#'
#' To be added
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references
#' \url{https://github.com/plbaldoni/epigraHMM}
#' 
#' @import Rhdf5lib
#' 
#' @examples 
#' 
#' # Creating dummy object
#' countData <- rbind(matrix(rnbinom(1e3,mu = 2,size = 10),ncol = 1),
#'                    matrix(rnbinom(2e3,mu = 7.5,size = 5),ncol = 1),
#'                    matrix(rnbinom(1e3,mu = 2,size = 10),ncol = 1))
#' 
#' colData <- data.frame(condition = 'A', replicate = 1)
#' object <- epigraHMMDataSetFromMatrix(countData,colData)
#' 
#' # Initializing
#' object <- initializer(object,controlEM())
#' 
#' # Running epigraHMM
#' object <- epigraHMM(object,controlEM(),type = 'consensus',dist = 'nb')
#'
#' @export
epigraHMM = function(object,control,type,dist = 'nb'){

    epsilon.em = criterion = NULL

    # Checking if the object is sorted
    if(!all(base::order(colData(object)[,c('condition','replicate')],decreasing = FALSE)==seq_len(nrow(colData(object))))){
        stop("Columns of epigraHMMDataSet are not sorted")
    }

    # Checking if the object is a RangedSummarizedExperiment
    if(!methods::is(object)[1]%in%c('SummarizedExperiment','RangedSummarizedExperiment')){
        stop("object is neither a SummarizedExperiment nor a RangedSummarizedExperiment")
    }

    # Check if control is correct
    if(!(all(names(control)%in%names(controlEM()))&all(names(controlEM())%in%names(control)))){
        stop("control is not a proper object. Use controlEM().")
    }

    # Check type
    if(!(type%in%c('consensus','differential'))){
        stop("type must be either 'consensus' or 'differential'")
    }

    # Check dist
    if(!(dist%in%c('zinb','nb'))){
        stop("dist must be either 'zinb' or 'nb'")
    }

    # Checking if type is appropriate
    if(type == 'consensus' & length(unique(colData(object)$condition))>1){
        stop("You requested type = 'consensus', but object has samples from multiple conditions. The input object must have data from a single condition in consensus peak calling.")
    }

    # Checking if type is appropriate
    if(type == 'differential' & length(unique(colData(object)$condition)) == 1){
        stop("You requested type = 'differential', but object has samples from a single condition. The input object must have data from multiple conditions in differential peak calling.")
    }

    # Running peak callers
    if(type == 'consensus'){
        object <- consensusHMM(object = object,control = control,dist = dist)
    } else{
        object <- differentialHMM(object = object,control = control,dist = dist)
    }
    return(object)
}
