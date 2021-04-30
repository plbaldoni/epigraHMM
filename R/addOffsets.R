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
#' @examples 
#' 
#' # Creating dummy object
#' countData <- list('counts' = matrix(rpois(4e5,10),ncol = 4),
#' 'controls' = matrix(rpois(4e5,5),ncol = 4))
#' colData <- data.frame(condition = c('A','A','B','B'), replicate = c(1,2,1,2))
#' object <- epigraHMMDataSetFromMatrix(countData,colData)
#'
#' # Adding pre-computed offsets
#' object <- addOffsets(object = object,
#'                      offsets = matrix(rnorm(4e5),ncol = 4))
#'
#' @export
addOffsets <- function(object, offsets) {
    
    if(any(is.nan(as.matrix(offsets)) | is.infinite(as.matrix(offsets)) | is.na(as.matrix(offsets)))){
        stop('offsets must not contain NA, NaN, or infinite values')
    }
    
    if ('offsets' %in% SummarizedExperiment::assayNames(object)) {
        SummarizedExperiment::assay(object, 'offsets') <-
            SummarizedExperiment::assay(object, 'offsets') + offsets
    } else{
        dimnames(offsets) <-
            dimnames(SummarizedExperiment::assay(object, 'counts'))
        SummarizedExperiment::assay(object, 'offsets', withDimnames = TRUE) <-
            offsets
    }
    return(object)
}