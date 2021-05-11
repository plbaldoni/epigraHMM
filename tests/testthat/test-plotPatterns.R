test_that("check if plot is working as it should",{
    
    # Creating dummy object
    countData <- cbind(rbind(matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
                             matrix(rnbinom(1e2, mu = 10, size = 5), ncol = 1),
                             matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
                             matrix(rnbinom(1e2, mu = 10, size = 5), ncol = 1),
                             matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
                             matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
                             matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1)),
                       rbind(matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
                             matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
                             matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
                             matrix(rnbinom(1e2, mu = 10, size = 5), ncol = 1),
                             matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1),
                             matrix(rnbinom(1e2, mu = 10, size = 5), ncol = 1),
                             matrix(rnbinom(1e2, mu = 1, size = 10), ncol = 1)))
    
    colData <- data.frame(condition = c('A','B'), replicate = c(1,1))
    rowRanges <- GenomicRanges::GRanges('chrA',
                                        IRanges::IRanges(start = seq(1,by = 500,
                                                                     length.out = nrow(countData)),width = 500))
    
    object <- epigraHMMDataSetFromMatrix(countData,colData,rowRanges = rowRanges)
    
    # Initializing
    object <- initializer(object,controlEM())
    
    # Running epigraHMM
    object <- epigraHMM(object,controlEM(maxIterEM = 2),type = 'differential',dist = 'nb')
    
    # Calling peaks
    peaks <- callPeaks(object = object,
                       hdf5 = S4Vectors::metadata(object)$output,
                       method = 'viterbi')
    
    # Plotting patterns
    fig <- plotPatterns(object,
                        ranges = peaks[1],
                        peaks = peaks)
    
    expect_equal(methods::is(fig),'pheatmap')
})
