################################################################################
################################################################################
### Functions w/ wrappers to get & generate objects directly from other packages
################################################################################
################################################################################

################################################################################
### Get fragment length
################################################################################

getFragLen <- function(bamFiles,gr.gaps,gr.blackList){
    gr.discard <- IRanges::union(gr.gaps, gr.blackList)
    myDiscard <- csaw::readParam(discard = gr.discard)
    fragList <- lapply(seq_len(length(bamFiles[['counts']])),function(x){
        fragCorr <- csaw::correlateReads(bam.files = bamFiles[['counts']][x], 
                                         param = myDiscard)
        return(csaw::maximizeCcf(fragCorr))
    })
    return(unlist(fragList))
}

################################################################################
### Get genome
################################################################################

getGenome <- function(genome,windowSize,bamFiles){
    if (methods::is(genome)[1] == "GRanges") {
        gr.genome <- unlist(GenomicRanges::tile(genome, width = windowSize))
        return(list('genome' = gr.genome, 'seqinfo' = NULL))
    } else{
        if (is.character(genome)) {
            chrlist <- lapply(bamFiles[['counts']],FUN = function(x){
                myPar <- Rsamtools::ScanBamParam(what = "rname")
                chrVal <- Rsamtools::scanBam(x, param = myPar)[[1]]$rname
                return(unique(as.character(chrVal)))
            })
            chrlist <- Reduce(intersect,chrlist)
            # Tiling up the specified genome
            gr.seqinfo <- GenomeInfoDb::Seqinfo(genome = genome)
            gr.genome <- GenomicRanges::tileGenome(cut.last.tile.in.chrom = TRUE,
                                                   seqlengths = gr.seqinfo,
                                                   tilewidth = windowSize)
            myDecode <- S4Vectors::decode(GenomeInfoDb::seqnames(gr.genome))
            gr.genome <- gr.genome[myDecode %in% chrlist]
            return(list('genome' = gr.genome,'seqinfo' = gr.seqinfo))
        } else{
            stop('The argument genome must be either a single string with the name of the reference genome (e.g. "hg19") or a "GRanges" object')
        }
    }
}

################################################################################
### Get blacklist
################################################################################

#' @importFrom GenomeInfoDb genome
getList <- function(blackList,genome){
    if (!methods::is(blackList)[1] == "GRanges") {
        greyItems <- utils::data(package = 'GreyListChIP')$results[,'Item']
        
        if (is.character(genome)) {
            dataExists <- paste0(genome,'.blacklist') %in% greyItems
        } else{
            genomeName <- unique(GenomeInfoDb::genome(genome))
            if (length(genomeName) == 1) {
                dataExists <- paste0(genomeName,'.blacklist') %in% greyItems
            } else{
                dataExists <- FALSE
            }
        }

        if (isTRUE(blackList) & is.character(genome) & dataExists) {
            greyName <- system.file(file.path('data',paste0(genome,'.blacklist.RData')),
                                    package = "GreyListChIP")
            load(greyName,envir = environment())
            gr.blacklist <- get(paste0(genome,'.blacklist'))
            gr.blackList <- 
                GenomicRanges::GRanges(seqnames = seqnames(gr.blacklist),
                                       ranges = IRanges::ranges(gr.blacklist),
                                       seqinfo = GenomeInfoDb::Seqinfo(genome = genome))
            gr.blackList <- GenomicRanges::trim(gr.blackList)
        } else {
            gr.blackList <- GenomicRanges::GRanges()
        }
    } else{
        gr.blackList <- blackList 
    }
    return(gr.blackList)
}

################################################################################
### Get UCSC gap track
################################################################################

getGap <- function(gapTrack,genome,gr.seqinfo = NULL){
    chrom <- NULL
    if (!methods::is(gapTrack)[1] == "GRanges") {
        if (isTRUE(gapTrack) & is.character(genome)) {
            # Gap table
            session <- rtracklayer::browserSession()
            GenomeInfoDb::genome(session) <- genome
            tb.ucsc <- rtracklayer::ucscTableQuery(session, table = "gap")
            dt.gaps <- rtracklayer::getTable(tb.ucsc)
            dt.gaps <- data.table::as.data.table(dt.gaps)
            gr.gaps <- GenomicRanges::makeGRangesFromDataFrame(
                df = dt.gaps[chrom %in% unique(seqnames(gr.seqinfo)), ],
                seqinfo = gr.seqinfo,
                starts.in.df.are.0based = TRUE)
        } else{
            gr.gaps <- GenomicRanges::GRanges()
        }
    } else{
        gr.gaps <- gapTrack
    }
    return(gr.gaps)
}

################################################################################
### Plot posterior probabilities
################################################################################

posteriorProbPlot <- function(fig.counts,DT,ifDifferential,hdf5,DTmelt,subsetIdx) {
    Sample = Window = P = name = ranges = .SD = NULL
    myColors <- c('Peaks' = '#4B9CD3','Annotation' = '#151515')
    myLabels <- c('Peaks' = ifelse(ifDifferential,'Differential Peaks','Consensus Peaks'),
                  'Annotation' = 'Annotation')
    xBreaks <- round(quantile(DT$Window, c(0.25, 0.5, 0.75), names = FALSE))
    xLabels <- scales::comma(DT$start[round(quantile(DT$Window, c(0.25, 0.5, 0.75),
                                                     names = FALSE))])
    if (!is.null(hdf5)) {
        dtProb <- cbind(DTmelt[Sample == levels(Sample)[1], c(.SD, 'name' = 'Prob.')],
                        'P' = exp(rhdf5::h5read(hdf5, 'logProb1')[subsetIdx, 2]))
        fig.counts <- fig.counts + 
            ggplot2::theme(axis.title.x = ggplot2::element_blank(),
                           axis.line.x = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank(),
                           axis.text.x = ggplot2::element_blank()) +
            ggplot2::scale_fill_manual(values = myColors,labels = myLabels) +
            ggplot2::scale_color_manual(values = myColors) +
            ggplot2::guides(color = 'none')
        fig.p <- ggplot2::ggplot(data = dtProb,ggplot2::aes(x = Window, y = P)) +
            ggplot2::geom_area(position = 'identity', alpha = 0.75,
                               color = '#4B9CD3',fill = '#4B9CD3') + 
            ggplot2::facet_grid(rows = ggplot2::vars(name)) +
            ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
            ggplot2::theme_bw() +
            ggplot2::labs(y = 'Prob.') +
            ggplot2::theme(strip.text.y = ggplot2::element_text(colour = ggplot2::alpha('grey', 0.0)),
                           legend.position = "none",panel.grid = ggplot2::element_blank()) +
            ggplot2::scale_x_continuous(breaks = xBreaks,
                                        labels = xLabels)
        if (!methods::is(ranges)[1] == "GRanges") {
            fig.p <- fig.p + ggplot2::labs(x = 'Genomic Window')
        } else {
            fig.p <- fig.p + ggplot2::labs(x = paste0('Genomic Window (', seqnames(ranges), ')'))
        }
        # Return
        fig <- ggpubr::ggarrange(fig.counts,fig.p,ncol = 1,nrow = 2,
                                 heights = c(0.8, 0.2),
                                 common.legend = TRUE,legend = 'top')
    } else{
        fig <- fig.counts + ggplot2::scale_x_continuous(breaks = xBreaks,labels = xLabels) +
            ggplot2::scale_fill_manual(values = myColors,labels = myLabels) +
            ggplot2::scale_color_manual(values = myColors) +
            ggplot2::guides(color = 'none')
        # Return
        if (!methods::is(ranges)[1] == "GRanges") {
            fig <- fig + ggplot2::labs(x = 'Genomic Window')
        } else {
            fig <- fig + ggplot2::labs(x = paste0('Genomic Window (', seqnames(ranges), ')'))
        }
    }
    return(fig)
}

################################################################################
### Plot counts
################################################################################

readCountPlot <- function(DTmelt,peaks,annotation){
    Window = Counts = Sample = .SD = NULL
    # Plotting counts
    fig.counts <- ggplot2::ggplot(data = DTmelt, ggplot2::aes(x = Window, y = Counts)) +
        ggplot2::facet_grid(rows = ggplot2::vars(Sample)) +
        ggplot2::geom_line() + ggplot2::theme_bw() +
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
    return(fig.counts)
}