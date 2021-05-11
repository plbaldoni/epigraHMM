#' Create a plot of differerential patterns posterior probabilities from epigraHMM
#'
#' `plotPatterns()` plots the posterior probabilities associated with differential patterns from a differential analysis of `epigraHMM()`
#'
#' @param object an epigraHMMDataSet
#' @param ranges a GRanges object or a pair of integers with the genomic corrdinates/windows to be plotted
#' @param peaks either a GRanges object or a vector of logicals (with length equal to the number of rows in `object`) specifying the genomic corrdinates/windows with peaks
#' @param hdf5 a character string with the hdf5 file path from `epigraHMM`
#' @param colors an optional argument that specifies the colors for each differential combinatorial pattern
#'
#' @return A pheatmat
#' 
#' @examples 
#' # Creating dummy object
#' countData <- cbind(rbind(matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 10, size = 5), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 10, size = 5), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1)),
#'                    rbind(matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 10, size = 5), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 10, size = 5), ncol = 1),
#'                          matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1)))
#' 
#' colData <- data.frame(condition = c('A','B'), replicate = c(1,1))
#' rowRanges <- GenomicRanges::GRanges('chrA',
#'                      IRanges::IRanges(start = seq(1,by = 500,
#'                      length.out = nrow(countData)),width = 500))
#' 
#' object <- epigraHMMDataSetFromMatrix(countData,colData,rowRanges = rowRanges)
#'
#' # Initializing
#' object <- initializer(object,controlEM())
#'
#' # Running epigraHMM
#' object <- epigraHMM(object,controlEM(),type = 'differential',dist = 'nb')
#'
#' # Calling peaks
#' peaks <- callPeaks(object = object,
#'                   hdf5 = S4Vectors::metadata(object)$output,
#'                   method = 'viterbi')
#' 
#' # Plotting patterns
#' plotPatterns(object,
#'              ranges = peaks[1],
#'              peaks = peaks)
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/epigraHMM}
#'
#' @importFrom pheatmap pheatmap
#' @importFrom grDevices colorRampPalette palette.colors
#'
#' @export
plotPatterns = function(object,
                        ranges,
                        peaks,
                        hdf5 = metadata(object)$output,
                        colors = NULL) {

    # Checking if differential analysis
    if (!(length(unique(object$condition)) > 1)) stop('Only one condition is present in the object input.')
    
    # Checking if hdf5 file exists
    if(!file.exists(hdf5)) stop('The hdf5 file does not exist.')
    
    # Reading in posterior probabilties and patterns
    postprob <- rhdf5::h5read(hdf5,'mixtureProb')
    patterns <- rhdf5::h5read(hdf5,'mixturePatterns')
    patterns <- unlist(lapply(patterns,function(x){paste(unique(object$condition)[as.numeric(gregexpr('E',x)[[1]])],collapse = '-')}))

    # Subsetting
    sub_index <- if (methods::is(ranges)[1] == "GRanges") overlapsAny(object, ranges) else seq(ranges[1], ranges[2])
    sub_object <- object[sub_index,]
    sub_postprob <- postprob[sub_index,]
    colnames(sub_postprob) <- patterns
    
    # Checking overlaps with differential peaks
    diff_index <- overlapsAny(sub_object, peaks)
    sub_postprob[!diff_index,] <- NA
    
    # Annotation 
    if(is.null(colors)) colors <- grDevices::palette.colors(n = 10, "Tableau 10")
    
    if(length(colors)<ncol(sub_postprob)) stop('Insufficient number of colors. Please provide ',ncol(sub_postprob),' different colors.')
        
    anno_colors <- list('Pattern' = colors[seq_len(ncol(sub_postprob))])
    names(anno_colors[['Pattern']]) <- colnames(sub_postprob)
    
    # Plotting
    fig <- pheatmap::pheatmap(sub_postprob,
                       color = grDevices::colorRampPalette(c(4, "white", 2))(256),
                       annotation_col = data.frame('Pattern' = colnames(sub_postprob),
                                                   row.names = colnames(sub_postprob)),
                       annotation_colors = anno_colors,
                       cluster_rows = FALSE,cluster_cols = FALSE,angle_col = 45,
                       main = paste0('Differential Enrichment Pattern\nPosterior Probability\n(',
                                     unique(seqnames(sub_object)),':',
                                     scales::comma(min(start(sub_object))),'-',
                                     scales::comma(max(end(sub_object))),')'),
                       labels_row = 'Genomic\nWindows',annotation_names_col = FALSE)
    
    return(fig)
}
