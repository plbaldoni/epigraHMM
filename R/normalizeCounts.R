#' Normalize counts
#'
#' This function performs a non-linear normalization of counts with respect to a reference sample (geometric mean)
#'
#' @param object an epigraHMMDataSet
#' @param control list of control arguments from controlEM()
#' @param ... arguments to be passed to \code{\link[limma]{loessFit}} for loess calculation
#'
#' @details
#'
#' This function `limma::loessFit`, which simply a wrapper for the `stats::lowess` smoother.
#'
#' @return An epigraHMMDataSet with an 'offsets' assay filled in.
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#'
#' @references
#' \url{https://github.com/plbaldoni/epigraHMM}
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment assayNames assay
#' @importFrom limma loessFit
#'
#' @export
normalizeCounts <- function(object,control,...){

    # Checking input
    if (!(methods::is(object)[1]=='RangedSummarizedExperiment')){
        stop('Check argments')
    }

    # Constructing reference sample
    if('offsets' %in% SummarizedExperiment::assayNames(object)){
        logSample <- log1p(SummarizedExperiment::assay(object,'counts'))-SummarizedExperiment::assay(object,'offsets')
    } else{
        logSample <- log1p(SummarizedExperiment::assay(object,'counts'))
    }
    logReference <- rowMeans(logSample)

    # Creating offsets
    offsets <- do.call(cbind,lapply(seq_len(ncol(SummarizedExperiment::assay(object,'counts'))),function(x){
        return(limma::loessFit(y=(logSample[,x] - logReference),x=(logSample[,x] + logReference)/2,...)$fitted)
    }))
    
    # Trimming offsets
    if(!is.null(control[['trimOffset']])){offsets <- round(offsets,digits = control[['trimOffset']])}

    # Organizing output
    if('offsets' %in% SummarizedExperiment::assayNames(object)){
        SummarizedExperiment::assay(object,'offsets') <- SummarizedExperiment::assay(object,'offsets') + offsets
    } else{
        dimnames(offsets) <- dimnames(SummarizedExperiment::assay(object,'counts'))
        SummarizedExperiment::assay(object,'offsets',withDimnames=TRUE) <- offsets
    }

    return(object)
}
