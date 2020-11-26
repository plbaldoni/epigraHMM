#' Create a epigraHMMDataSet from a set of BAM files
#'
#' This function creates a \code{\link[SummarizedExperiment]{RangedSummarizedExperiment}} object from of a set of BAM files.
#' It is used to store the input data, the model offsets, and the results from the peak calling algorithms.
#'
#' @param bamFiles a string vector (or a list of string vectors) with the path for BAM files. If bamFiles is a list of string vectors,
#' vectors must be named, have the same dimension, and, at least, a vector with name 'counts' must exist (see details).
#' @param colData a \code{data.frame} with the experimental data. It must contain the columns \code{condition} and \code{replicate}.
#' \code{condition} refers to the experimental condition identifier (e.g. cell line name). \code{replicate} refers to the replicate identification number (unique for each condition).
#' @param genome either a single string with the name of the reference genome (e.g. 'hg19') or a GRanges object with ranges to be tilled into a set of non-overlapping windows.
#' @param windowSize an integer specifying the size of genomic windows where read counts will be computed.
#' @param gapTrack either a logical (\code{TRUE}, the default, or \code{FALSE}) or a GRanges object with gap regions of the genome to be excluded. If \code{TRUE}, the function will discard genomic coordinates overlapping regions present in the UCSC gap table of the respective reference genome (if available). See Details section below.
#' @param blackList either a logical (\code{TRUE}, the default, or \code{FALSE}) or a GRanges object with blacklisted regions of the genome to be excluded. If \code{TRUE}, the function will discard ENCODE blacklisted regions from selected reference genomes (if available). See Details section below.
#'
#' @details
#'
#' The index ".bai" files must be stored in the same directory of their respective BAM files.
#' The index files must be named after their respective BAM files with the additional ".bai" suffix.
#' 
#' `epigraHMMDataSetFromBam` will store experimental data (e.g. ChIP-seq counts) from bamFiles (or bamFiles[['counts']], if a list is provided).
#' Additional data (e.g. input control counts) will be stored similarly with their respective list names.
#' 
#' By default, the function computes read counts using csaw's estimated fragment length via cross correlation analysis.
#' For experimental counts (e.g. ChIP-seq), sequencing reads are shifted downstream half of the estimated fragment length.
#' For additional counts (e.g. input control), sequencing reads are not shifted prior to counting.
#' 
#' Additional columns included in the colData input will be passed to the resulting epigraHMMDataSet assay and can be acessed via \code{colData()} function.
#'
#' The \code{genome} argument will call GenomeInfoDb::Seqinfo() to fetch the chromosome lengths of the specified genome.
#' See ?GenomeInfoDb::Seqinfo for the list of UCSC genomes that are currently supported.
#'
#' If \code{gapTrack = TRUE} and the name of a reference genome is passed as input through \code{genome} (e.g. 'hg19'),
#' the function will discard any genomic coordinate overlapping regions specified by the respective UCSC gap table.
#' If \code{gapTrack} is a GRanges object, the function will discard any genomic coordinate overlaping regions from \code{gapTrack}.
#'
#' If \code{blackList = TRUE} and the name of a reference genome is passed as input through \code{genome} (e.g. 'hg19'),
#' The function will fetch the manually curated blacklist tracks (Version 2) from \url{https://github.com/Boyle-Lab/Blacklist/tree/master/lists}.
#' Current available genomes are ce10, dm3, hg19, hg38, and mm10.
#' If \code{blackList} is a GRanges object, the function will discard any genomic coordinate overlaping regions from \code{blackList}.
#'
#' @return An epigraHMMDataSet object with sorted colData regarding conditions and replicates.
#' Experimental counts will be stored in the 'counts' assay in the resulting epigraHMMDataSet object.
#' Additional experimental data will be stored with their respective names from the list bamFiles.
#'
#' @author Pedro L. Baldoni, \email{pedrobaldoni@gmail.com}
#'
#' @references
#' \url{https://github.com/plbaldoni/epigraHMM}
#' DOI: 10.1093/nar/gkv1191
#' DOI: 10.1038/s41598-019-45839-z
#' DOI: 10.1038/nature11247
#'
#' @importFrom methods is
#' @importFrom SummarizedExperiment SummarizedExperiment rowRanges assay colData
#' @importFrom GenomeInfoDb Seqinfo seqnames genome
#' @importFrom GenomicRanges tile tileGenome makeGRangesFromDataFrame GRanges
#' @importFrom rtracklayer browserSession getTable ucscTableQuery import BEDFile
#' @importFrom IRanges overlapsAny union
#' @importFrom Rsamtools scanBam ScanBamParam
#' @importFrom bamsignals bamCount
#' @importFrom csaw maximizeCcf correlateReads readParam
#' @importFrom data.table as.data.table
#' @importFrom S4Vectors decode
#' 
#' @examples
#' \dontrun{
#' Rsamtools::indexBam(MMDiffBamSubset::WT.AB2())
#' Rsamtools::indexBam(MMDiffBamSubset::Input())
#' 
#' bamFiles <- list('counts' = MMDiffBamSubset::WT.AB2(),
#'                 'controls' = MMDiffBamSubset::Input())
#'                 
#' colData <- data.frame(condition = 'WT',replicate = 1)
#'                 
#' epigraHMMDataSetFromBam(bamFiles,colData,genome = 'mm10',windowSize = 250,
#'                         gapTrack = TRUE,
#'                         blackList = TRUE)
#' }
#'
#' @export
epigraHMMDataSetFromBam <- function(bamFiles,
                                    colData,
                                    genome,
                                    windowSize,
                                    gapTrack = TRUE,
                                    blackList = TRUE){
    
    condition = replicate = chrom = NULL
    
    # Checking if colData has the correct format
    if(!(is.data.frame(colData) & all(c('condition','replicate')%in%names(colData)))){
        stop("The argument colData must be a data.frame with the columns 'condition' and 'replicate'")
    }
    
    # Checking whether replicates are unique
    if(any(table(colData$condition)>1)){
        if(nrow(unique(colData[,c('condition','replicate')]))<nrow(colData[,c('condition','replicate')])){
            stop('The columns "condition" and "replicate" must uniquely represent your data')
        }
    }
    
    # Checking whether bamFiles is correct
    if(is.list(bamFiles) & ('counts' %in% names(bamFiles))){
        if(any(unlist(lapply(bamFiles,function(x){
            return(!(is.character(x) & length(x)==nrow(colData) & all(file.exists(x)) & all(file.exists(paste0(x,'.bai')))))
        })))){
            stop("bamFiles is not a proper argument, check the help manual.")
        }
    } else{
        if(!(is.character(bamFiles) & length(bamFiles)==nrow(colData) & all(file.exists(bamFiles)) & all(file.exists(paste0(bamFiles,'.bai'))))){
            stop("bamFiles is not a proper argument, check the help manual.")
        } else{
            names(bamFiles) <- 'counts'
        }
    }
    
    # Checking whether windowSize is an integer
    if(!(is.numeric(windowSize) & windowSize%%1==0)){
        stop('The argument windowSize must be an integer number')
    }
    
    # Setting up reference genome
    if(methods::is(genome)[1]=="GRanges"){
        gr.genome <- unlist(GenomicRanges::tile(genome,width = windowSize))
    } else{
        if(is.character(genome)){
            chrlist <- Reduce(intersect,lapply(bamFiles[['counts']],FUN = function(x){unique(as.character(Rsamtools::scanBam(x,param = Rsamtools::ScanBamParam(what = "rname"))[[1]]$rname))}))

            # Tiling up the specified genome
            gr.seqinfo <- GenomeInfoDb::Seqinfo(genome=genome)
            gr.genome <- GenomicRanges::tileGenome(seqlengths = gr.seqinfo,tilewidth = windowSize,cut.last.tile.in.chrom = TRUE)
            gr.genome <- gr.genome[S4Vectors::decode(GenomeInfoDb::seqnames(gr.genome))%in%chrlist]
        } else{
            stop('The argument genome must be either a single string with the name of the reference genome (e.g. "hg19") or a "GRanges" object')
        }
    }
    
    # Setting up gap track
    if(gapTrack == TRUE & is.character(genome)){
        # Gap table
        session <- rtracklayer::browserSession()
        GenomeInfoDb::genome(session) <- genome
        gr.gaps <- GenomicRanges::makeGRangesFromDataFrame(
            df = data.table::as.data.table(rtracklayer::getTable(rtracklayer::ucscTableQuery(session,table="gap")))[chrom%in%unique(seqnames(gr.seqinfo)),],
            seqinfo = gr.seqinfo,
            starts.in.df.are.0based = TRUE)
    } else{
        if(methods::is(gapTrack)[1]=="GRanges"){
            gr.gaps <- gapTrack
        } else{
            gr.gaps <- GenomicRanges::GRanges()
        }
    }
    
    # Setting up blacklist track
    if(blackList == TRUE & is.character(genome)){
        gr.blackList <- GenomicRanges::trim(GenomicRanges::makeGRangesFromDataFrame(df = data.table::as.data.table(rtracklayer::import(rtracklayer::BEDFile(paste0('https://github.com/Boyle-Lab/Blacklist/raw/master/lists/',genome,'-blacklist.v2.bed.gz')))),
                                                                                    seqinfo = GenomeInfoDb::Seqinfo(genome=genome)))
    } else{
        if(methods::is(blackList)[1]=="GRanges"){
            gr.blackList <- blackList
        } else{
            gr.blackList <- GenomicRanges::GRanges()
        }
    }
    
    # Cleaning up the genome
    gr.genome <- gr.genome[!IRanges::overlapsAny(gr.genome,IRanges::union(gr.gaps,gr.blackList))]
    
    # Estimating the fragment length
    colData$fragLength <- unlist(lapply(seq_len(nrow(colData)),function(x){
        csaw::maximizeCcf(csaw::correlateReads(bam.files = bamFiles[['counts']][x],param=csaw::readParam(discard=IRanges::union(gr.gaps,gr.blackList))))
    }))
    
    # Computing read counts and adding to the output
    epigraHMMDataSet <- SummarizedExperiment::SummarizedExperiment(assays = list(counts = matrix(do.call(cbind,lapply(seq_len(nrow(colData)),FUN = function(x){return(bamsignals::bamCount(bampath = bamFiles[['counts']][x],gr = gr.genome,verbose = FALSE,shift = colData[['fragLength']][x]/2))})),
                                                                                                 byrow = FALSE,nrow = length(gr.genome),ncol = nrow(colData),dimnames = list(NULL,paste(colData$condition,colData$replicate,sep='.')))),
                                                                   rowRanges = gr.genome,
                                                                   colData = colData)
    
    
    # If there are controls, repeat
    if(is.list(bamFiles)){
        for(idx in names(bamFiles)[-which(names(bamFiles)=='counts')]){
            tmp <- do.call(cbind,lapply(seq_len(nrow(colData)),FUN = function(x){return(bamsignals::bamCount(bampath = bamFiles[[idx]][x],gr = gr.genome,verbose = FALSE))}))
            dimnames(tmp) <- dimnames(SummarizedExperiment::assay(epigraHMMDataSet,'counts'))
            SummarizedExperiment::assay(epigraHMMDataSet,idx) <- tmp
        }
        rm(tmp)
    }
    
    # Sorting the object
    if(!all(base::order(SummarizedExperiment::colData(epigraHMMDataSet)[,c('condition','replicate')],decreasing = F)==seq_len(nrow(SummarizedExperiment::colData(epigraHMMDataSet))))){
        epigraHMMDataSet <- epigraHMMDataSet[,base::order(SummarizedExperiment::colData(epigraHMMDataSet)[,c('condition','replicate')],decreasing = F)]
        message("Rows of colData have been sorted with respect to conditions (and replicates). The resulting colData is:")
        message(paste0(utils::capture.output(SummarizedExperiment::colData(epigraHMMDataSet)), collapse = "\n"))
    }
    
    return(epigraHMMDataSet)
}