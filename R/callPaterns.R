#' Extract posterior probabilities (or combinatorial patterns) associated with differential regions
#'
#' Given results from epigraHMM's differential peak caller, this function will 
#' output either posterior probabilities or combinatorial patterns associated 
#' with the mixture components of the embedded mixture model.
#' 
#' @param object an epigraHMMDataSet
#' @param peaks a GRanges object with differential peaks from `callPeaks`
#' @param hdf5 a character with the location of the epigraHMM HDF5 output file
#' @param type a character string that defines which output will be givem 
#' (see details; default is 'all')
#' @param fdr the desired fdr thresholding level to define combinatorial 
#' patterns
#' @param pattern a string that explicitly specifies the combinatorial pattern
#' to be output
#' @param ranges a GRanges object with the genomic ranges to subset the output
#' 
#' @details 
#' 
#' The output of `callPatterns` is always restricted to genomic windows
#' intersecting peaks.
#' 
#' If `type = 'all'`, all windows' posterior probabilities associated with
#' all differential combinatorial patterns are returned. If `type = 'fdr'`,
#' users must also specify the input argument `pattern` and this function will
#' output windows wich are associated with the given `pattern` that pass a
#' particular fdr threshold level. If `type = 'max'`, this function will output
#' the combinatorial pattern which has the maximal posterior probability for 
#' each window. If `type = 'ranges'`, the windows that are output are restricted
#' to those that intersect the `ranges` input argument.
#' 
#' @return A GRanges object with metadata
#' 
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references \url{https://github.com/plbaldoni/epigraHMM}
#' 
#' @examples 
#' # Creating dummy object
#' countData <- cbind(rbind(matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 10, size = 5), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 10, size = 5), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1)),
#'                    rbind(matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 10, size = 5), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 10, size = 5), ncol = 1),
#'                          matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1)))
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
#' # Extracting posterior probabilities
#' patterns <- callPatterns(object = object,peaks = peaks,type = 'max')
#'
#' @importFrom S4Vectors DataFrame mcols
#' @export
callPatterns <- function(object,peaks,hdf5 = metadata(object)$output,
                         type = 'all',fdr = NULL,pattern = NULL,ranges = NULL){
    
    out <- rowRanges(object)
    
    postprob <- data.table(rhdf5::h5read(hdf5,'mixtureProb'))
    colnames(postprob) <- getPatterns(x = object, hdf5 =  hdf5)
    
    idx <- overlapsAny(out,peaks)
    out <- out[idx]
    postprob <- postprob[idx,]
    
    if(type == 'all'){
        S4Vectors::mcols(out) <- DataFrame(postprob)
    }
    if(type == 'fdr'){
        ppFdr <- DataFrame(fdrControl(prob = postprob[[pattern]],fdr = fdr))
        colnames(ppFdr) <- pattern
        S4Vectors::mcols(out) <- DataFrame(ppFdr)
    }
    if(type == 'max'){
        ppmax <- postprob[,list(Enrichment =colnames(postprob)[which.max(.SD)]),
                          by = seq_len(nrow(postprob))]
        
        S4Vectors::mcols(out) <- DataFrame(ppmax[,'Enrichment'])
    }
    if(type == 'ranges'){
        idxRanges <- overlapsAny(out,ranges)
        out <- out[idxRanges]
        S4Vectors::mcols(out) <- DataFrame(postprob[idxRanges,])
    }
    
    return(out)
}