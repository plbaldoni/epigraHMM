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
#' @export
initializer <- function(object,control){

    # Checking input
    if (!(methods::is(object)[1]=='RangedSummarizedExperiment')){
        stop('Check argments')
    }
    if('viterbi' %in% SummarizedExperiment::assayNames(object)){
        stop('viterbi assay already exists')
    }

    # Organizing output
    viterbi <- Matrix::Matrix(do.call(cbind,lapply(seq_len(ncol(assay(object))),function(x){
        initializerHMM(object = object[,x],control = control)
    })),dimnames = dimnames(assay(object)),sparse = TRUE)

    SummarizedExperiment::assay(object,'peaks',withDimnames=TRUE) <- viterbi

    return(object)
}
