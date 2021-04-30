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
    if (methods::is(ranges)[1] == "GRanges") {
        subsetIdx <- overlapsAny(object, ranges)
    } else{
        subsetIdx <- seq(ranges[1], ranges[2])
    }
    subobject <- object[subsetIdx,]
    
    mat <- as.matrix(SummarizedExperiment::assay(subobject, 'counts') / exp(SummarizedExperiment::assay(subobject, 'offsets')))
    DT <- data.table::as.data.table(mat)
    rm(mat)
    
    # Transforming the data (consensus or differential?)
    ifDifferential <- (length(unique(object$condition)) > 1)
    
    DT <- aggregateCounts(DT,ifDifferential,object)
    
    # Melting data.table
    list2env(meltCounts(DT,ranges,object,subobject,peaks,annotation,subsetIdx),
             envir = environment())
    
    # Plotting counts
    fig.counts <- ggplot2::ggplot(data = DTmelt, ggplot2::aes(x = Window, y = Counts)) +
        ggplot2::facet_grid(rows = ggplot2::vars(Sample)) +
        ggplot2::geom_line() +
        ggplot2::theme_bw() +
        ggplot2::labs(y = 'Normalized Read Counts', x = 'Genomic Window') +
        ggplot2::theme(panel.grid = ggplot2::element_blank())
    
    # Adding peak track
    if (!is.null(peaks)) {
        DT_anno <- DTmelt[Sample == levels(Sample)[1], c(.SD, 'name' = 'Peaks')]
        fig.counts <- fig.counts + 
            ggplot2::geom_rect(data = DT_anno,ggplot2::aes(xmin = Window,
                                                           xmax = Window,
                                                           ymin = 0.99 * (max(DTmelt$Counts) + 5) * ifelse(peaks, peaks, NA),
                                                           ymax = 1.01 * (max(DTmelt$Counts) + 5) * ifelse(peaks, peaks, NA),
                                                           color = 'Peaks',fill = 'Peaks'),na.rm = TRUE) +
            ggplot2::theme(legend.position = 'top',
                           legend.title = ggplot2::element_blank(),
                           legend.direction = 'horizontal')
    }
    
    # Adding annotation
    if (!is.null(annotation)) {
        DT_anno <- DTmelt[Sample == levels(Sample)[1], c(.SD, 'name' = 'Peaks')]
        fig.counts <- fig.counts + 
            ggplot2::geom_rect(data = DT_anno,ggplot2::aes(xmin = Window,
                                                           xmax = Window,
                                                           ymin = 0.99 * (max(DTmelt$Counts) + 10) * ifelse(annotation, annotation, NA),
                                                           ymax = 1.01 * (max(DTmelt$Counts) + 10) * ifelse(annotation, annotation, NA),
                                                           color = 'Annotation',fill = 'Annotation'),na.rm = TRUE) +
            ggplot2::theme(legend.position = 'top',
                           legend.title = ggplot2::element_blank(),
                           legend.direction = 'horizontal')
    }
    
    # Adding plotting posterior probabilities
    if (!is.null(hdf5)) {
        
        dtProb <- cbind(DTmelt[Sample == levels(Sample)[1], c(.SD, 'name' = 'Prob.')],
                        'P' = exp(rhdf5::h5read(hdf5, 'logProb1')[subsetIdx, 2]))
        
        fig.counts <- fig.counts + 
            ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                           axis.line.x = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank()) +
            ggplot2::scale_fill_manual(values = c('Peaks' = '#4B9CD3','Annotation' = '#151515'),
                                       labels = c('Peaks' = ifelse(ifDifferential,'Differential Peaks','Consensus Peaks'),'Annotation' = 'Annotation')) +
            ggplot2::scale_color_manual(values = c('Peaks' = '#4B9CD3','Annotation' = '#151515')) +
            ggplot2::guides(color = FALSE)
        
        fig.p <- ggplot2::ggplot(data = dtProb,ggplot2::aes(x = Window, y = P)) +
            ggplot2::geom_area(position = 'identity', alpha = 0.75,
                               color = '#4B9CD3',fill = '#4B9CD3') + 
            ggplot2::facet_grid(rows = ggplot2::vars(name)) +
            ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
            ggplot2::theme_bw() +
            ggplot2::labs(y = 'Prob.') +
            ggplot2::theme(strip.text.y = ggplot2::element_text(colour = ggplot2::alpha('grey', 0.0)),
                           legend.position = "none",
                           panel.grid = ggplot2::element_blank()) +
            ggplot2::scale_x_continuous(breaks = round(quantile(DT$Window, c(0.25, 0.5, 0.75), names = FALSE)),
                                        labels = scales::comma(DT$start[round(quantile(DT$Window, c(0.25, 0.5, 0.75), names = FALSE))]))
        
        if (!methods::is(ranges)[1] == "GRanges") {
            fig.p <- fig.p + 
                ggplot2::labs(x = 'Genomic Window')
        } else{
            fig.p <- fig.p + 
                ggplot2::labs(x = paste0('Genomic Window (', seqnames(ranges), ')'))
        }
        
        # Return
        fig <- ggpubr::ggarrange(fig.counts,
                                 fig.p,
                                 ncol = 1,
                                 nrow = 2,
                                 heights = c(0.8, 0.2),
                                 common.legend = TRUE,
                                 legend = 'top')
    } else{
        # Return
        fig <- fig.counts + ggplot2::scale_x_continuous(breaks = round(quantile(DT$Window, c(0.25, 0.5, 0.75), names = FALSE)),
                                                        labels = scales::comma(DT$start[round(quantile(DT$Window, c(0.25, 0.5, 0.75), names = FALSE))])) +
            ggplot2::scale_fill_manual(values = c('Peaks' = '#4B9CD3','Annotation' = '#151515'),
                                       labels = c('Peaks' = ifelse(ifDifferential,'Differential Peaks','Consensus Peaks'),'Annotation' = 'Annotation')) +
            ggplot2::scale_color_manual(values = c('Peaks' = '#4B9CD3','Annotation' = '#151515')) +
            ggplot2::guides(color = FALSE)
        
        if (!methods::is(ranges)[1] == "GRanges") {
            fig <- fig +
                ggplot2::labs(x = 'Genomic Window')
        } else{
            fig <- fig +
                ggplot2::labs(x = paste0('Genomic Window (', seqnames(ranges), ')'))
        }
    }
    
    return(fig)
}
