test_that("check viterbi path from initialization",{
    
    # Creating dummy object
    set.seed(210423)
    countData <- rbind(matrix(rnbinom(2e3,mu = 2,size = 10),ncol = 2),
                       matrix(rnbinom(4e3,mu = 100,size = 5),ncol = 2),
                       matrix(rnbinom(2e3,mu = 2,size = 10),ncol = 2))
    
    colData <- data.frame(condition = c('A','A'), replicate = c(1,2))
    
    rowRanges <- GenomicRanges::GRanges('chrA',
                                        IRanges::IRanges(start = seq(from = 1, length.out = 4e3,by = 250),width = 250))
    
    object <- epigraHMMDataSetFromMatrix(countData,colData,rowRanges)
    
    # Initializing
    object <- initializer(object,controlEM())
    
    expect_true('peaks' %in% SummarizedExperiment::assayNames(object))
    expect_equal(dim(assay(object,'peaks')),dim(countData))
    expect_equal(unique(c(as.matrix(assay(object,'peaks')))),c(0,1))
})
