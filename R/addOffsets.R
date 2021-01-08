#' Add offsets to epigraHMMDataSet
#'
#' This function adds model offsets to epigraHMMDataSet
#'
#' @param object an epigraHMMDataSet
#' @param offsets a matrix with model offsets
#'
#' @details
#'
#' To be added
#'
#' @return An epigraHMMDataSet with an 'offsets' assay filled in.
#'
#' @references
#' \url{https://github.com/plbaldoni/epigraHMM}
#'
#' @importFrom SummarizedExperiment assayNames assay
#'
#' @export
addOffsets <- function(object,offsets){
    if('offsets' %in% SummarizedExperiment::assayNames(object)){
        SummarizedExperiment::assay(object,'offsets') <- SummarizedExperiment::assay(object,'offsets') + offsets
    } else{
        dimnames(offsets) <- dimnames(SummarizedExperiment::assay(object,'counts'))
        SummarizedExperiment::assay(object,'offsets',withDimnames=TRUE) <- offsets
    }
    return(object)
}