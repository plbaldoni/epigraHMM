################################################################################
################################################################################
### R script w/ functions to write the output of epigraHMM to file
################################################################################
################################################################################

################################################################################
### Function to save output files
################################################################################

#' @importFrom S4Vectors subjectHits
saveOutputFiles <- function(gr.bed,object,control,hdf5,method) {
    subjectHits = peakindex = NULL
    # Adding necessary columns for browser
    prob <- exp(h5read(hdf5,'logProb1')[,2])
    gr.overlaps <- SummarizedExperiment::findOverlaps(gr.bed,rowRanges(object))
    dt.overlaps <- data.table::as.data.table(gr.overlaps)
    dt.overlaps <- dt.overlaps[,mean(prob[subjectHits]),by = 'queryHits']
    gr.bed$score <- 1000*dt.overlaps$V1
    gr.bed$thickStart <- start(gr.bed)
    gr.bed$thickEnd <- end(gr.bed)
    # File names
    chrset <- as.character(unique(seqnames(rowRanges(object))))
    dirPath <- path.expand(control[['tempDir']])
    fileName <- paste0(control[['fileName']],'_',
                       c('peaks.bed',paste0('prob_',chrset,'.wig')))
    filePath <- file.path(dirPath,fileName)
    filenames <- vapply(filePath,checkPath,FUN.VALUE = 'character',
                        USE.NAMES = FALSE)
    names(filenames) <- c('peaks',paste0('prob_',chrset))
    # Writing bed file with peaks
    dt.bed <- data.frame(chrom = seqnames(gr.bed), chromStart = start(gr.bed), 
                         chromEnd = end(gr.bed), name = gr.bed$name, 
                         score = gr.bed$score, strand = '.', 
                         thickStart = gr.bed$thickStart, 
                         thickEnd = gr.bed$thickEnd)
    writeBed(dt.bed,control,method,filenames)
    # Writing wig files with posterior probabilities
    dt.bigwig <- data.frame(chrom = seqnames(rowRanges(object)),
                            chromStart = start(rowRanges(object)),
                            chromEnd = end(rowRanges(object)), prob = prob)
    writeWig(object,chrset,dt.bigwig,control,filenames)
    rm(dt.bigwig)
    # Writing bedGraph files for mixture probabilities
    if (length(unique((colData(object)[['condition']]))) > 1) {
        mixProbSet <- h5read(hdf5,'mixturePatterns')
        fileName <- paste0(control[['fileName']],'_',
                           paste0('mixProb_',mixProbSet,'.wig'))
        wigPath <- file.path(dirPath,fileName)
        filenames <- 
            c(filenames,vapply(wigPath,checkPath,FUN.VALUE = 'character',
                               USE.NAMES = FALSE))
        names(filenames)[names(filenames) == ""] <- mixProbSet
        for (i in seq_len(length(mixProbSet))) {
            mixProb <- h5read(hdf5,'mixtureProb')[peakindex,i]
            dt.bedgraph = data.frame(chrom = seqnames(gr.graph),
                                     chromStart = start(gr.graph),
                                     chromEnd = end(gr.graph),
                                     mixProb = mixProb)
            writeBedGraph(dt.bed,dt.bedgraph,control,mixProbSet[i],filenames)
        }
    }
    message('The following files have been saved:')
    for (i in seq_len(length(filenames))) {
        message(filenames[i])
    }
}

################################################################################
### Write bedGraph file
################################################################################

writeBedGraph <- function(dt.bed,dt.bedgraph,control,mixProbSet,filenames){
    outPath <- control[['tempDir']]
    utils::write.table(dt.bedgraph,
                       file = file.path(outPath, "temp.bed"), row.names = FALSE,
                       col.names = FALSE, quote = FALSE, sep = "\t")
    header1 = paste0('track type=bedGraph name="epigraHMM(mixProb:',mixProbSet,')" description="','Probability','" visibility=full maxHeightPixels=128:32:11 graphType=bar autoScale=off alwaysZero=on viewLimits=0.0:1.0')
    header2 = paste0('browser position ',dt.bed[1,'chrom'],':',
                     dt.bed[1,'chromStart'],'-',dt.bed[1,'chromEnd'])
    system2('echo',paste0(header2,' | cat - ',file.path(outPath,"temp.bed"),
                          ' > ',file.path(outPath,"temp1.bed")))
    system2('echo',paste0(header1,' | cat - ',file.path(outPath,"temp1.bed"),
                          ' > ',as.character(filenames[mixProbSet])))
    system2('rm',paste(file.path(outPath,"temp.bed"),
                       file.path(outPath,"temp1.bed")))
}

################################################################################
### Write WIG file
################################################################################

writeWig <- function(object,chrset,dt.bigwig,control,filenames){
    outPath <- control[['tempDir']]
    for (i in chrset) {
        dt.bigwig.subset <- dt.bigwig[dt.bigwig$chrom == i, ]
        
        utils::write.table(dt.bigwig.subset[,c('chromStart','prob')],
                           file = file.path(outPath, "temp.bed"),
                           row.names = FALSE,col.names = FALSE, quote = FALSE, 
                           sep = "\t")
        header1 = paste0('track type=wiggle_0 name="epigraHMM(prob:',i,')" description="','Probability','" visibility=full maxHeightPixels=128:32:11 graphType=bar autoScale=off alwaysZero=on viewLimits=0.0:1.0')
        header2 = paste0('browser position ',dt.bigwig.subset[1,'chrom'],':',
                         dt.bigwig.subset[1,'chromStart'],'-',
                         dt.bigwig.subset[10,'chromStart'])
        header3 = paste0('variableStep chrom=',i)
        system2('echo',paste0(header3,' | cat - ',file.path(outPath,"temp.bed"),
                              ' > ',file.path(outPath,"temp1.bed")))
        system2('echo',paste0(header1,' | cat - ',file.path(outPath,"temp1.bed"),
                              ' > ',file.path(outPath,"temp2.bed")))
        system2('echo',paste0(header2,' | cat - ',file.path(outPath,"temp2.bed"),
                              ' > ',as.character(filenames[paste0('prob_',i)])))
        system2('rm',paste(file.path(outPath,"temp.bed"),
                           file.path(outPath,"temp1.bed"),
                           file.path(outPath,"temp2.bed")))
        tryCatch({
            sqInfo <- 
                GenomeInfoDb::seqinfo(SummarizedExperiment::rowRanges(object))
            rtracklayer::wigToBigWig(x = filenames[paste0('prob_',i)],
                                     seqinfo = sqInfo)
            system2("rm",filenames[paste0('prob_',i)])
        },error = function(x){
            message("It was not possible to convert wig files to BigWig format because the input object has no specified genome")
        })
    }
}

################################################################################
### Write BED file
################################################################################

writeBed <- function(dt.bed,control,method,filenames){
    outPath <- control[['tempDir']]
    utils::write.table(dt.bed, file = file.path(outPath, "temp.bed"),
                       row.names = FALSE, col.names = FALSE, quote = FALSE,
                       sep = "\t")
    header1 <- paste0('track name="epigraHMM" description="',ifelse(method == 'viterbi','Viterbi Peaks',paste0('FDR-controlled Peaks (FDR = ',method,')')),'" visibility=1 useScore=1')
    header2 <- paste0('browser position ',dt.bed[1,'chrom'],':',
                      dt.bed[1,'chromStart'],'-',dt.bed[1,'chromEnd'])
    system2('echo',paste0(header2,' | cat - ',file.path(outPath,"temp.bed"),
                          ' > ',file.path(outPath,"temp1.bed")))
    system2('echo',paste0(header1,' | cat - ',file.path(outPath,"temp1.bed"),
                          ' > ',as.character(filenames['peaks'])))
    system2('rm',paste(file.path(outPath,"temp.bed"),
                       file.path(outPath,"temp1.bed")))
}