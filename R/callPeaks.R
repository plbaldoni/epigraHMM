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
callPeaks = function(object,
                     hdf5,
                     method = 'viterbi',
                     saveToFile = FALSE,
                     control = NULL)
{
    subjectHits = i = NULL

    # Calling peaks
    prob <- exp(rhdf5::h5read(hdf5,'logProb1')[,2])
    if(method=='viterbi'){
        peakindex <- (rhdf5::h5read(hdf5,'viterbi')[,1]==1)
    } else{
        if(is.numeric(method) & method>0 & method<1){
            peakindex <- fdrControl(prob = prob,fdr = method)
        } else{
            stop('The argument method is not valid')
        }
    }
    gr.graph <- SummarizedExperiment::rowRanges(object)[peakindex]
    gr.bed <- GenomicRanges::reduce(gr.graph)

    # Summarize the output
    gr.bed$name <- paste0(paste0('peak',seq_len(length(gr.bed))))
    gr.bed$score <- 1000*data.table::as.data.table(IRanges::findOverlaps(gr.bed,SummarizedExperiment::rowRanges(object)))[,mean(prob[subjectHits]),by='queryHits']$V1
    gr.bed$thickStart <- GenomicRanges::start(gr.bed)
    gr.bed$thickEnd <- GenomicRanges::end(gr.bed)

    # File names
    if(saveToFile){
        chrset <- as.character(unique(SummarizedExperiment::seqnames(SummarizedExperiment::rowRanges(object))))

        filenames <- sapply(file.path(path.expand(control[['tempDir']]),paste0(control[['fileName']],'_',c('peaks.bed',paste0('prob_',chrset,'.wig')))),checkPath,USE.NAMES = FALSE)
        names(filenames) <- c('peaks',paste0('prob_',chrset))

        # Writing bed file with peaks
        dt.bed <- data.frame(chrom=SummarizedExperiment::seqnames(gr.bed),chromStart=SummarizedExperiment::start(gr.bed),
                             chromEnd=SummarizedExperiment::end(gr.bed),name=gr.bed$name,score=gr.bed$score,strand='.',
                             thickStart=gr.bed$thickStart,thickEnd=gr.bed$thickEnd)

        utils::write.table(dt.bed,file=file.path(control[['tempDir']],"temp.bed"),row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")
        header1 <- paste0('track name="epigraHMM" description="',ifelse(method == 'viterbi','Viterbi Peaks',paste0('FDR-controlled Peaks (FDR = ',method,')')),'" visibility=1 useScore=1')
        header2 <- paste0('browser position ',dt.bed[1,'chrom'],':',dt.bed[1,'chromStart'],'-',dt.bed[1,'chromEnd'])

        system2('echo',paste0(header2,' | cat - ',file.path(control[['tempDir']],"temp.bed"),' > ',file.path(control[['tempDir']],"temp1.bed")))
        system2('echo',paste0(header1,' | cat - ',file.path(control[['tempDir']],"temp1.bed"),' > ',as.character(filenames['peaks'])))
        system2('rm',paste(file.path(control[['tempDir']],"temp.bed"),file.path(control[['tempDir']],"temp1.bed")))

        # Writing wig files with posterior probabilities
        dt.bigwig <- data.frame(chrom=SummarizedExperiment::seqnames(SummarizedExperiment::rowRanges(object)),
                                chromStart=SummarizedExperiment::start(SummarizedExperiment::rowRanges(object)),
                                chromEnd=SummarizedExperiment::end(SummarizedExperiment::rowRanges(object)),
                                prob=prob)

        for(i in chrset){
            dt.bigwig.subset <- dt.bigwig[dt.bigwig$chrom==i,]

            utils::write.table(dt.bigwig.subset[,c('chromStart','prob')],
                               file=file.path(control[['tempDir']],"temp.bed"),row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")
            header1 = paste0('track type=wiggle_0 name="epigraHMM(prob:',i,')" description="','Probability','" visibility=full maxHeightPixels=128:32:11 graphType=bar autoScale=off alwaysZero=on viewLimits=0.0:1.0')
            header2 = paste0('browser position ',dt.bigwig.subset[1,'chrom'],':',dt.bigwig.subset[1,'chromStart'],'-',dt.bigwig.subset[10,'chromStart'])
            header3 = paste0('variableStep chrom=',i)

            system2('echo',paste0(header3,' | cat - ',file.path(control[['tempDir']],"temp.bed"),' > ',file.path(control[['tempDir']],"temp1.bed")))
            system2('echo',paste0(header1,' | cat - ',file.path(control[['tempDir']],"temp1.bed"),' > ',file.path(control[['tempDir']],"temp2.bed")))
            system2('echo',paste0(header2,' | cat - ',file.path(control[['tempDir']],"temp2.bed"),' > ',as.character(filenames[paste0('prob_',i)])))
            system2('rm',paste(file.path(control[['tempDir']],"temp.bed"),file.path(control[['tempDir']],"temp1.bed"),file.path(control[['tempDir']],"temp2.bed")))

            tryCatch({
                rtracklayer::wigToBigWig(x = filenames[paste0('prob_',i)],seqinfo = GenomeInfoDb::seqinfo(SummarizedExperiment::rowRanges(object)))
                system2("rm",filenames[paste0('prob_',i)])
            },error = function(x){message("It was not possible to convert wig files to BigWig format because the input object has no specified genome")})
        }
        rm(dt.bigwig)

        # Writing bedGraph files for mixture probabilities
        if(length(unique((colData(object)[['condition']])))>1){

            mixProbSet <- rhdf5::h5read(hdf5,'mixturePatterns')
            filenames <- c(filenames,sapply(file.path(path.expand(control[['tempDir']]),paste0(control[['fileName']],'_',paste0('mixProb_',mixProbSet,'.wig'))),checkPath,USE.NAMES = FALSE))
            names(filenames)[names(filenames)==""] <- mixProbSet

            for(i in seq_len(length(mixProbSet))){

                dt.bedgraph = data.frame(chrom=SummarizedExperiment::seqnames(gr.graph),chromStart=SummarizedExperiment::start(gr.graph),
                                         chromEnd=SummarizedExperiment::end(gr.graph),
                                         mixProb=rhdf5::h5read(hdf5,'mixtureProb')[peakindex,i])

                utils::write.table(dt.bedgraph,file=file.path(control[['tempDir']],"temp.bed"),row.names = FALSE,col.names = FALSE,quote = FALSE,sep="\t")
                header1 = paste0('track type=bedGraph name="epigraHMM(mixProb:',mixProbSet[i],')" description="','Probability','" visibility=full maxHeightPixels=128:32:11 graphType=bar autoScale=off alwaysZero=on viewLimits=0.0:1.0')
                header2 = paste0('browser position ',dt.bed[1,'chrom'],':',dt.bed[1,'chromStart'],'-',dt.bed[1,'chromEnd'])

                system2('echo',paste0(header2,' | cat - ',file.path(control[['tempDir']],"temp.bed"),' > ',file.path(control[['tempDir']],"temp1.bed")))
                system2('echo',paste0(header1,' | cat - ',file.path(control[['tempDir']],"temp1.bed"),' > ',as.character(filenames[mixProbSet[i]])))
                system2('rm',paste(file.path(control[['tempDir']],"temp.bed"),file.path(control[['tempDir']],"temp1.bed")))
            }
        }

        message('The following files have been saved:')
        for(i in seq_len(length(filenames))){message(filenames[i])}
    }

    return(gr.bed)
}
