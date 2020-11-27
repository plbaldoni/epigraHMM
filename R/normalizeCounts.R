#' Normalize counts
#'
#' This function performs a non-linear normalization of counts with respect to a reference sample (geometric mean)
#'
#' @param object an epigraHMMDataSet
#' @param type either 'loess', the default, or 'gam'
#' @param ... arguments to be passed to \code{\link[mgcv]{gam}} or \code{\link[limma]{loessFit}} for loess calculation
#'
#' @details
#'
#' If `type = 'gam'`, `normalizeCounts` will use `mgcv::bam` with option `discrete = TRUE`.
#' If `type = 'loess'`, it will use `limma::loessFit`, which simply a wrapper for the `stats::lowess` smoother.
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
#' @importFrom mgcv bam s
#'
#' @export
normalizeCounts <- function(object,type = 'loess',...){

    #@importFrom mgcv bam s predict.bam

    # Checking input
    if (!(methods::is(object)[1]=='RangedSummarizedExperiment') & (type %in% c('gam','loess'))){
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
        if (type == 'loess') {
            return(limma::loessFit(y=(logSample[,x] - logReference),x=(logSample[,x] + logReference)/2,...)$fitted)
        } else{
            return(mgcv::bam(y~s(x,bs = 'cs'),discrete = TRUE,data = data.frame(y = (logSample[,x] - logReference),
                                                                                x = (logSample[,x] + logReference)/2))$fitted.values)
        }
    }))

    # Organizing output
    if('offsets' %in% SummarizedExperiment::assayNames(object)){
        SummarizedExperiment::assay(object,'offsets') <- SummarizedExperiment::assay(object,'offsets') + offsets
    } else{
        dimnames(offsets) <- dimnames(SummarizedExperiment::assay(object,'counts'))
        SummarizedExperiment::assay(object,'offsets',withDimnames=TRUE) <- offsets
    }

    return(object)
}
