#' Create a plot with the results from epigraHMM
#'
#' `plotCounts()` plots read counts and peak regions from `epigraHMM()`
#'
#' @param object an epigraHMMDataSet
#' @param ranges a GRanges object or a pair of integers with the genomic corrdinates/windows to be plotted
#' @param peaks an optional parameter with a GRanges object or a vector of logicals (with length equal to the number of rows in `object`) specifying the genomic corrdinates/windows with peaks
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
                     peaks = NULL,
                     hdf5 = NULL) {
    Sample = Counts = name = P = NULL

    # Checking ranges
    if (!((methods::is(ranges)[1] == "GRanges" &
           length(ranges) == 1) |
          (methods::is(ranges)[1] == "numeric" &
           length(ranges) == 2)
    )) {
        stop(
            'The argument ranges must be either a GRanges object with length one or a numeric vector of integers of length two'
        )
    }

    # Subset input
    if (methods::is(ranges)[1] == "GRanges") {
        subsetIdx <- overlapsAny(object, ranges)
    } else{
        subsetIdx <- seq(ranges[1], ranges[2])
    }
    subobject <- object[subsetIdx, ]


    DT <-
        data.table::as.data.table(as.matrix(
            SummarizedExperiment::assay(subobject, 'counts') / exp(SummarizedExperiment::assay(subobject, 'offsets'))
        ))

    # Transforming the data (consensus or differential?)
    ifDifferential <-
        (length(unique(
            SummarizedExperiment::colData(object)$condition
        )) > 1)
    if (ifDifferential) {
        nameCol <-
            paste0(unique(SummarizedExperiment::colData(object)$condition))
        for (i in nameCol) {
            DT[, paste0(i) := rowSums(.SD), .SDcols = which(SummarizedExperiment::colData(object)$condition ==
                                                                i)]
        }
        DT <- DT[, nameCol, with = FALSE]
    } else{
        nameCol <-
            paste0('Replicate ',
                   SummarizedExperiment::colData(object)$replicate)
        data.table::setnames(DT, nameCol)
    }

    # Bringing in extra variables
    if (methods::is(ranges)[1] == "GRanges") {
        DT <-
            cbind(DT, as.data.table(rowRanges(object))[subsetIdx, c('seqnames', 'start', 'end', 'width', 'strand')])
        DTmelt <-
            data.table::melt(
                DT,
                id.vars = c('seqnames', 'start', 'end'),
                measure.vars = seq_len(ncol(DT) - 5),
                value.name = 'Counts',
                variable.name = 'Sample'
            )
        if (!is.null(peaks)) {
            DTmelt[Sample == levels(Sample)[1], peaks := overlapsAny(subobject, peaks)]
        }
    } else{
        DT[, start := seq(ranges[1], ranges[2])]
        DTmelt <-
            data.table::melt(
                DT,
                id.vars = c('start'),
                measure.vars = seq_len(ncol(DT) - 1),
                value.name = 'Counts',
                variable.name = 'Sample'
            )
        if (!is.null(peaks)) {
            DTmelt[Sample == levels(Sample)[1], peaks := peaks[subsetIdx]]
        }
    }

    # Plotting counts
    fig.counts <-
        ggplot2::ggplot(data = DTmelt, ggplot2::aes(x = start, y = Counts)) +
        ggplot2::facet_grid(rows = ggplot2::vars(Sample)) +
        ggplot2::geom_line() +
        ggplot2::theme_bw() + ggplot2::labs(y = 'Normalized Read Counts', x = 'Genomic Window') +
        ggplot2::theme(panel.grid = ggplot2::element_blank())

    if (!is.null(peaks)) {
        fig.counts <-
            fig.counts + ggplot2::geom_path(
                data = DTmelt[Sample == levels(Sample)[1], c(.SD, 'name' = 'Peaks')],
                lineend = "square",
                ggplot2::aes(
                    x = start,
                    y = (max(DTmelt$Counts) + 5) * ifelse(peaks, peaks, NA),
                    color = name,
                ),
                na.rm = TRUE,
                size = 5,
            ) +
            ggplot2::scale_color_manual(
                values = c('#4B9CD3'),
                labels = ifelse(
                    ifDifferential,
                    'Differential Peaks',
                    'Consensus Peaks'
                )
            ) +
            ggplot2::theme(
                legend.position = 'top',
                legend.title = ggplot2::element_blank(),
                legend.direction = 'horizontal'
            ) +
            ggplot2::scale_x_continuous(labels = scales::comma)
    }

    # Plotting posterior probabilities
    if (!is.null(hdf5)) {
        fig.counts <- fig.counts +
            ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.line.x = ggplot2::element_blank(),
                axis.ticks.x = ggplot2::element_blank(),
                axis.text.x = ggplot2::element_blank()
            )


        fig.p <-
            ggplot2::ggplot(data = cbind(DTmelt[Sample == levels(Sample)[1], c(.SD, 'name' = 'Prob.')], 'P' = exp(
                rhdf5::h5read(hdf5, 'logProb1')[subsetIdx, 2]
            )),
            ggplot2::aes(x = start, y = P)) +
            ggplot2::facet_grid(rows = ggplot2::vars(name)) +
            ggplot2::geom_area(
                position = 'identity',
                alpha = 0.75,
                color = '#4B9CD3',
                fill = '#4B9CD3'
            ) +
            ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
            ggplot2::scale_x_continuous(labels = scales::comma) +
            ggplot2::theme_bw() +
            ggplot2::labs(y = 'Prob.') +
            ggplot2::theme(
                strip.text.y = ggplot2::element_text(colour = ggplot2::alpha('grey', 0.0)),
                legend.position = "none",
                panel.grid = ggplot2::element_blank()
            )

        if (!methods::is(ranges)[1] == "GRanges") {
            fig.p <- fig.p + ggplot2::labs(x = 'Genomic Window')
        } else{
            fig.p <-
                fig.p + ggplot2::labs(x = paste0('Genomic Window (', seqnames(ranges), ')'))
        }

        # Return
        fig <-
            ggpubr::ggarrange(
                fig.counts,
                fig.p,
                ncol = 1,
                nrow = 2,
                heights = c(0.8, 0.2),
                common.legend = TRUE,
                legend = 'top'
            )
    } else{
        fig <- fig.counts
    }

    return(fig)
    # # Now, plotting
    # if (plot == 'counts') {
    #     # Peak track
    #     maxcount <-
    #         max(datamelt[Sample == unique(datamelt$Sample)[1] &
    #                          plotindex == TRUE, Counts])
    #     dataanno = data.table::data.table(
    #         as.data.table(rowRanges(object)),
    #         Sample = unique(datamelt$Sample)[1],
    #         plotindex = plotindex
    #     )
    #     if (type == 'viterbi') {
    #         dataanno[, Peak := (S4Vectors::metadata(object)$viterbi == ifelse(is.consensus, 'E', 'D'))]
    #     } else{
    #         if (is.numeric(type) & type > 0 & type < 1) {
    #             dataanno[, Peak := fdrControl(S4Vectors::metadata(object)$prob[[colindex]], fdr = type)]
    #         } else{
    #             dataanno[, Peak := FALSE]
    #         }
    #     }
    #     dataanno$Peak = plyr::mapvalues(dataanno$Peak,
    #                                     from = c(FALSE, TRUE),
    #                                     to = c(NA, colindex))
    #     dataanno[, PeakCount := ifelse(is.na(Peak), NA, 1.1 * maxcount)]
    #     dataanno[, PeakTrack := ifelse(colindex == 'Differential',
    #                                    'Differential Peak Call',
    #                                    'Consensus Peak Call')]
    #
    #     # Gene track
    #     if (length(unique(genome(object))) == 1 & annotate) {
    #         # Selecting genes
    #         session <- rtracklayer::browserSession()
    #         GenomeInfoDb::genome(session) <-
    #             unique(genome(object))[1]
    #         grGene <- GenomicRanges::makeGRangesFromDataFrame(
    #             rtracklayer::getTable(
    #                 rtracklayer::ucscTableQuery(session, table = "ncbiRefSeq")
    #             ),
    #             starts.in.df.are.0based = TRUE,
    #             seqnames.field = 'chrom',
    #             start.field = 'txStart',
    #             end.field = 'txEnd',
    #             strand.field = 'strand',
    #             keep.extra.columns = TRUE
    #         )
    #
    #         # Reducing each gene based on their transcripts start/end sites
    #         grGene <- unlist(reduce(split(grGene, grGene$name2)))
    #
    #         # Overlapping genes with rowRanges
    #         dataanno[, Gene := overlapsAny(rowRanges(object), grGene)]
    #
    #         dataanno$Gene = plyr::mapvalues(
    #             dataanno$Gene,
    #             from = c(FALSE, TRUE),
    #             to = c(NA, colindex)
    #         )
    #         dataanno[, GeneCount := ifelse(is.na(Gene), NA, 1.2 * maxcount)]
    #         dataanno[, GeneTrack := 'NCBI RefSeq Genes']
    #
    #         # Creating gene name track
    #         refseq.out <-
    #             as.data.table(grGene)[, c('seqnames', 'start', 'end')][, gene_name := names(grGene)]
    #         refseq.out[, plotindex := overlapsAny(grGene, with(unique(datamelt[(plotindex == TRUE), c('seqnames', 'start', 'end')]), GRanges(
    #             seqnames, IRanges(start, end)
    #         )))]
    #         refseq.out[plotindex == TRUE, Counts := 1.2 * maxcount]
    #         refseq.out[plotindex == TRUE, Lbl.Counts := 1.2 * maxcount]
    #         refseq.out[, Sample := unique(datamelt$Sample)[1]]
    #
    #         refseq.out[(plotindex == TRUE), start := max(c(start, start(ranges))), by = "gene_name"]
    #         refseq.out[(plotindex == TRUE), end := min(c(end, end(ranges))), by = "gene_name"]
    #
    #         refseq.out[plotindex == TRUE, Lbl.x := start + 0.5 * (end -
    #                                                                   start), by = "gene_name"]
    #     }
    #
    #     colval <- c('#4B9CD3', '#151515')
    #     names(colval) <-
    #         c(unique(dataanno$PeakTrack), 'NCBI RefSeq Genes')
    #
    #     fig.ChIP <-
    #         ggplot2::ggplot(data = datamelt[(plotindex == TRUE), ], ggplot2::aes(x =
    #                                                                                  start, y = Counts)) +
    #         ggplot2::geom_line() +
    #         ggplot2::facet_grid(rows = ggplot2::vars(Sample)) +
    #         ggplot2::geom_segment(
    #             inherit.aes = FALSE,
    #             data = dataanno[(plotindex == TRUE) &
    #                                 !is.na(Peak), ],
    #             ggplot2::aes(
    #                 x = start,
    #                 xend = end,
    #                 y = PeakCount,
    #                 yend = PeakCount,
    #                 color = PeakTrack
    #             ),
    #             size = 2
    #         ) +
    #         ggplot2::theme_bw() + ggplot2::labs(y = 'Normalized Read Counts') +
    #         ggplot2::scale_color_manual(values = colval) +
    #         ggplot2::theme(
    #             axis.title.x = ggplot2::element_blank(),
    #             axis.line.x = ggplot2::element_blank(),
    #             axis.ticks.x = ggplot2::element_blank(),
    #             axis.text.x = ggplot2::element_blank(),
    #             legend.position = 'top',
    #             legend.direction = 'horizontal',
    #             legend.title = ggplot2::element_blank(),
    #             panel.grid = ggplot2::element_blank()
    #         )
    #
    #     if (annotate) {
    #         fig.ChIP <- fig.ChIP +
    #             ggplot2::geom_segment(
    #                 inherit.aes = FALSE,
    #                 data = dataanno[(plotindex == TRUE) &
    #                                     !is.na(Gene), ],
    #                 ggplot2::aes(
    #                     x = start,
    #                     xend = end,
    #                     y = GeneCount,
    #                     yend = GeneCount,
    #                     color = GeneTrack
    #                 ),
    #                 size = 2
    #             ) +
    #             ggrepel::geom_text_repel(
    #                 data = refseq.out[(plotindex == TRUE), ],
    #                 ggplot2::aes(
    #                     x = Lbl.x,
    #                     y = Counts,
    #                     label = gene_name
    #                 ),
    #                 vjust = -5,
    #                 segment.size = 0.2,
    #                 size = 2,
    #                 segment.color = 'grey',
    #                 direction = 'y'
    #             ) +
    #             ggplot2::scale_y_continuous(limits = c(NA, 1.20 * max(dataanno$GeneCount, na.rm = T)))
    #     }
    #
    #     ### Figure 2: Post. Probabilities ###
    #     PostProb <-
    #         data.table::data.table(
    #             as.data.table(rowRanges(object)),
    #             Label = 'Prob',
    #             Prob = S4Vectors::metadata(object)$prob[[colindex]],
    #             plotindex = plotindex
    #         )
    #
    #     fig.Prob <-
    #         ggplot2::ggplot(data = PostProb[plotindex == TRUE, ], ggplot2::aes(x =
    #                                                                                start, y = Prob)) +
    #         ggplot2::facet_grid(rows = ggplot2::vars(Label)) +
    #         ggplot2::geom_area(
    #             position = 'identity',
    #             alpha = 0.5,
    #             color = '#4B9CD3',
    #             fill = '#4B9CD3'
    #         ) +
    #         ggplot2::scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
    #         ggplot2::scale_x_continuous(labels = scales::comma) +
    #         ggplot2::labs(x = paste0('Genomic Window (', seqnames(ranges), ')'),
    #                       y = 'Prob.') +
    #         ggplot2::theme_bw() +
    #         ggplot2::theme(strip.text.y = ggplot2::element_text(colour = ggplot2::alpha('grey', 0.0))) +
    #         ggplot2::theme(legend.position = "none",
    #                        panel.grid = ggplot2::element_blank())
    #
    #     fig <-
    #         ggpubr::ggarrange(
    #             fig.ChIP,
    #             fig.Prob,
    #             ncol = 1,
    #             nrow = 2,
    #             heights = c(0.8, 0.2),
    #             common.legend = TRUE,
    #             legend = 'top'
    #         )
    #     return(fig)
    # }
}
