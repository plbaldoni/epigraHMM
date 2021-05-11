#' Create a plot with the results from epigraHMM
#'
#' `plotCounts()` plots read counts and peak regions from `epigraHMM()`
#'
#' @param object an epigraHMMDataSet
#' @param ranges a GRanges object or a pair of integers with the genomic corrdinates/windows to be plotted
#' @param peaks an optional parameter with a GRanges object or a vector of logicals (with length equal to the number of rows in `object`) specifying the genomic corrdinates/windows with peaks
#' @param annotation an optional parameter with a GRanges object or a vector of logicals (with length equal to the number of rows in `object`) specifying the genomic corrdinates/windows of an annotation track
#' @param hdf5 an optional character string with the hdf5 file path from `epigraHMM`
#'
#' @details
#'
#' If the input object contains the assay 'offset',
#' reads will be normalized prior to plotting (e.g. counts/exp(offset)).
#' Reads from replicates pertaining to the same condition are aggregated prior to plotting.
#'
#' @return A ggplot
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/epigraHMM}
#'
#' @importFrom ggpubr ggarrange
#' @importFrom scales comma
#' @importFrom ggplot2 ggplot aes facet_grid geom_line theme_bw labs theme element_blank geom_path scale_color_manual scale_x_continuous geom_area scale_y_continuous element_text
#' @importFrom IRanges overlapsAny
#' @importFrom stats quantile
#'
#' @examples
#'
#' countData <- rbind(matrix(rnbinom(1e3,mu = 2,size = 10),ncol = 1),
#'                    matrix(rnbinom(1e3,mu = 7.5,size = 5),ncol = 1),
#'                    matrix(rnbinom(1e3,mu = 7.5,size = 5),ncol = 1),
#'                    matrix(rnbinom(1e3,mu = 2,size = 10),ncol = 1))
#'
#' colData <- data.frame(condition = 'A', replicate = 1)
#'
#' object <- epigraHMMDataSetFromMatrix(countData,colData)
#'
#' plotCounts(object,ranges = c(500,3500))
#'
#' @export
plotCounts = function(object,
                      ranges,
                      hdf5 = metadata(object)$output,
                      peaks = NULL,
                      annotation = NULL) {
    Sample = Counts = name = P = Window = DTmelt =  NULL
    
    # Checking ranges
    if (!((methods::is(ranges)[1] == "GRanges" & length(ranges) == 1) | 
          (methods::is(ranges)[1] %in% c("integer", "numeric") & length(ranges) == 2))) {
        stop('The argument ranges must be either a GRanges object with length one or a numeric vector of integers of length two')
    }
    
    # If peaks/annotation is not NULL, then make sure data types are consistent
    checkPlot(x = peaks, ranges = ranges)
    checkPlot(x = annotation, ranges = ranges)
    
    # Subset input
    subsetIdx <- if (methods::is(ranges)[1] == "GRanges") overlapsAny(object, ranges) else seq(ranges[1], ranges[2])
    subobject <- object[subsetIdx,]
    
    mat <- as.matrix(SummarizedExperiment::assay(subobject, 'counts') / exp(SummarizedExperiment::assay(subobject, 'offsets')))
    DT <- data.table::as.data.table(mat)
    rm(mat)
    
    # Transforming the data (consensus or differential?)
    ifDifferential <- (length(unique(object$condition)) > 1)
    
    DT <- aggregateCounts(DT,ifDifferential,object)
    
    # Melting data.table
    list2env(meltCounts(DT,ranges,object,subobject,peaks,annotation,subsetIdx),envir = environment())
    
    # Plotting counts
    fig.counts <- readCountPlot(DTmelt,peaks,annotation)
    
    # Adding plotting posterior probabilities
    fig <- posteriorProbPlot(fig.counts,DT,ifDifferential,hdf5,DTmelt,subsetIdx)
    
    return(fig)
}
