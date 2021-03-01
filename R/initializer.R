#' Initializer of epigraHMM
#'
#' This function call enriched windows individually for each sample in an epigraHMMDataSet.
#' These are then used for initializing purposes in epigraHMM. By default, the Viterbi algorithm
#' is used to determine enriched windows. Input controls and normalizing offsets are not utilized
#' in this initialization step.
#'
#' @param object an epigraHMMDataSet
#' @param control list of control arguments from controlEM()
#'
#' @details
#'
#' To be added
#'
#' @return An epigraHMMDataSet with a 'peaks' assay filled in.
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#'
#' @references
#' \url{https://github.com/plbaldoni/epigraHMM}
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom Matrix Matrix
#' 
#' @examples 
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
#' # Visualizing initialization peaks
#' #plot(assay(object),type = 'l')
#' #lines(7.5*assay(object,'peaks'),col = 'red')
#'
#' @export
initializer <- function(object,control){

    # Checking input
    if(!methods::is(object)[1]%in%c('SummarizedExperiment','RangedSummarizedExperiment')){
        stop("object is neither a SummarizedExperiment nor a RangedSummarizedExperiment")
    }
    if('peaks' %in% SummarizedExperiment::assayNames(object)){
        stop('initializing assay already exists')
    }

    # Organizing output
    viterbi <- Matrix::Matrix(do.call(cbind,lapply(seq_len(ncol(assay(object))),function(x){
        initializerHMM(object = object[,x],control = control)
    })),dimnames = dimnames(assay(object)),sparse = TRUE)

    SummarizedExperiment::assay(object,'peaks',withDimnames=TRUE) <- viterbi

    return(object)
}
