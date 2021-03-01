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
#' @importFrom data.table data.table
#' 
#' @examples
#' 
#' # Creating dummy object
#' countData <- list('counts' = matrix(rpois(1e5,10),ncol = 2),
#' 'controls' = matrix(rpois(1e5,5),ncol = 2))
#' colData <- data.frame(condition = c('A','A'), replicate = c(1,2))
#' object <- epigraHMMDataSetFromMatrix(countData,colData)
#'
#' # Normalizing counts
#' object <- normalizeCounts(object = object,control = controlEM(), span = 1)
#'
#' @export
normalizeCounts <- function(object,control,...){
    
    weights = minus = average = NULL

    # Checking input
    if(!methods::is(object)[1]%in%c('SummarizedExperiment','RangedSummarizedExperiment')){
        stop("object is neither a SummarizedExperiment nor a RangedSummarizedExperiment")
    }

    # Constructing reference sample
    if('offsets' %in% SummarizedExperiment::assayNames(object)){
        logSample <- log1p(SummarizedExperiment::assay(object,'counts'))-SummarizedExperiment::assay(object,'offsets')
    } else{
        logSample <- log1p(SummarizedExperiment::assay(object,'counts'))
    }
    logReference <- Matrix::rowMeans(logSample)

    # Creating offsets
    offsets <- do.call(cbind,lapply(seq_len(ncol(SummarizedExperiment::assay(object,'counts'))),function(x){
        
        dt_all <- data.table::data.table(minus = (logSample[,x] - logReference),
                             average = (logSample[,x] + logReference)/2,
                             weights = 1)
        
        dt_reduced <- dt_all[,list(weights = sum(weights)),by = c('minus','average')][,offsets := limma::loessFit(y=minus,x=average,weights=weights,min.weight = 0,max.weight = Inf,...)$fitted]
        
        dt_all <- merge(dt_all,dt_reduced[,c('minus','average','offsets')],by = c('minus','average'),all.x = TRUE, sort = FALSE)
        
        return(dt_all$offsets)
    }))
    
    # Trimming offsets
    if(!is.null(control[['trimOffset']])){offsets <- round(offsets,digits = control[['trimOffset']])}

    # Organizing output
    object <- addOffsets(object,offsets)

    return(object)
}
