#' Summarize peak calls and optionally create a BED 6+3 file in broadPeak format for visualization
#'
#' This function imports the output from `epigraHMM` and outputs a set of
#' peaks (consensus or differential) for a given FDR control threshold or Viterbi sequence.
#'
#' @param object an epigraHMMDataSet
#' @param hdf5 a character with the location of the epigraHMM HDF5 output file
#' @param method either 'viterbi' or a numeric FDR control threshold (e.g. 0.05). Default is 'viterbi'.
#' @param saveToFile a logical indicating whether or not to save the results to file.
#' Output files are always saved with peaks of interest defined on the region level. Default is FALSE.
#' @param control list of control arguments from controlEM(). This is an optional parameter and it is
#' only required when `saveToFile = TRUE` so that the output directory can be obtained. Default is NULL.
#'
#' @return A GRanges object with differential peak calls in BED 6+3 format
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/epigraHMM}
#'
#' @importFrom S4Vectors metadata
#' @importFrom rhdf5 h5read
#' @importFrom SummarizedExperiment rowRanges seqnames
#' @importFrom GenomicRanges reduce start end
#' @importFrom data.table as.data.table
#' @importFrom rtracklayer wigToBigWig
#' 
#' @examples 
#' 
#' # Creating dummy object
#' countData <- rbind(matrix(rnbinom(1e3,mu = 2,size = 10),ncol = 1),
#'                    matrix(rnbinom(2e3,mu = 7.5,size = 5),ncol = 1),
#'                    matrix(rnbinom(1e3,mu = 2,size = 10),ncol = 1))
#' 
#' colData <- data.frame(condition = 'A', replicate = 1)
#' 
#' rowRanges <- GenomicRanges::GRanges('chrA',
#' IRanges::IRanges(start = seq(from = 1, length.out = 4e3,by = 250),width = 250))
#' 
#' object <- epigraHMMDataSetFromMatrix(countData,colData,rowRanges)
#' 
#' # Initializing
#' object <- initializer(object,controlEM())
#' 
#' # Running epigraHMM
#' object <- epigraHMM(object,controlEM(),type = 'consensus',dist = 'nb')
#' 
#' # Calling peaks
#' peaks <- callPeaks(object = object,
#'                    hdf5 = S4Vectors::metadata(object)$output,
#'                    method = 'viterbi')
#'
#' @export
callPeaks = function(object,hdf5 = metadata(object)$output,method = 'viterbi',
                     saveToFile = FALSE,control = NULL)
{
    subjectHits = i = NULL
    
    # Checking if hdf5 exists
    if(!file.exists(hdf5)) stop('Input hdf5 file does not exist')

    # Calling peaks
    prob <- exp(h5read(hdf5,'logProb1')[,2])
    peakindex <- if(method=='viterbi') (h5read(hdf5,'viterbi')[,1]==1) else fdrControl(prob = prob,fdr = method)
    
    # If there is no rowRanges, return vector
    if(is.null(rowRanges(object))) return(peakindex)
    
    gr.graph <- rowRanges(object)[peakindex]
    gr.bed <- reduce(gr.graph)

    # Summarize the output
    gr.bed$name <- paste0(paste0('peak',seq_len(length(gr.bed))))

    # File names
    if(saveToFile) saveOutputFiles(gr.bed,object,control,hdf5)

    return(gr.bed)
}
