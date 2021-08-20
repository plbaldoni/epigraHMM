#' Segmentation of a genome in non-overlapping windows
#' 
#' This function segments a genome into non-overlapping windows.
#'
#' @param genome a string with the name of the genome (e.g. 'hg19')
#' @param window an integer with the window size
#' @param rm.gap a logical indicating gap regions should be removed
#' @param rm.blacklist a logical indicating blacklisted regions should be 
#' removed
#' 
#' @return a GRanges object with the binned genome
#' 
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#' @references
#' \url{https://github.com/plbaldoni/epigraHMM}
#' 
#' @examples
#' 
#' gr <- segmentGenome(genome = 'mm10', window = 500)
#'
#' @export
segmentGenome <- function(genome, window, rm.gap = TRUE, rm.blacklist = TRUE){
    
    # Getting genome info
    genomeSeqinfo <- Seqinfo(genome = genome)
    
    # Getting gap and blacklist track
    grGaps <- 
        getGap(gapTrack = rm.gap, genome = genome,gr.seqinfo = genomeSeqinfo)
    
    grBlacklist <- getList(blackList = rm.blacklist, genome = genome)
    
    # Getting genome bins
    gr.genome <- 
        tileGenome(cut.last.tile.in.chrom = TRUE,
                   tilewidth = window,seqlengths = genomeSeqinfo)
    
    #Excluding gap and blacklisted regions
    gr.genome <- gr.genome[!overlapsAny(gr.genome,union(grGaps,grBlacklist))]
    
    return(gr.genome)
}