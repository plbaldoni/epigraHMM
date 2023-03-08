test_that("check consensus (NB)",{
    
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
    
    # Running epigraHMM
    object <- epigraHMM(object,controlEM(),type = 'consensus',dist = 'nb')
    
    # Calling peaks with Viterbi
    peaks <- callPeaks(object = object,
                       hdf5 = S4Vectors::metadata(object)$output,
                       method = 'viterbi')
    
    # Checking boundaries
    tb <- table(call = overlapsAny(object,peaks),truth = c(rep(FALSE,1e3),rep(TRUE,2e3),rep(FALSE,1e3)))
    expect_true(tb[2,1]/sum(tb[2,])<=0.05)
    
    # Calling peaks with FDR thresholding
    peaks <- callPeaks(object = object,
                       hdf5 = S4Vectors::metadata(object)$output,
                       method = 0.05)
    
    # Checking boundaries
    tb <- table(call = overlapsAny(object,peaks),truth = c(rep(FALSE,1e3),rep(TRUE,2e3),rep(FALSE,1e3)))
    expect_true(tb[2,1]/sum(tb[2,])<=0.05)
})

test_that("check consensus (ZINB)",{
    
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
    
    # Running epigraHMM
    object <- epigraHMM(object,controlEM(),type = 'consensus',dist = 'zinb')
    
    # Calling peaks with Viterbi
    peaks <- callPeaks(object = object,
                       hdf5 = S4Vectors::metadata(object)$output,
                       method = 'viterbi')
    
    # Checking boundaries
    tb <- table(call = overlapsAny(object,peaks),truth = c(rep(FALSE,1e3),rep(TRUE,2e3),rep(FALSE,1e3)))
    expect_true(tb[2,1]/sum(tb[2,])<=0.05)
    
    # Calling peaks with FDR thresholding
    peaks <- callPeaks(object = object,
                       hdf5 = S4Vectors::metadata(object)$output,
                       method = 0.05)
    
    # Checking boundaries
    tb <- table(call = overlapsAny(object,peaks),truth = c(rep(FALSE,1e3),rep(TRUE,2e3),rep(FALSE,1e3)))
    expect_true(tb[2,1]/sum(tb[2,])<=0.05)
})

test_that("check output files from callPeaks",{
  
  # Creating dummy object
  set.seed(210423)
  countData <- cbind(rbind(matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
                           matrix(rnbinom(1e3, mu = 10, size = 5), ncol = 1),
                           matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
                           matrix(rnbinom(1e3, mu = 10, size = 5), ncol = 1),
                           matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
                           matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
                           matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1)),
                     rbind(matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
                           matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
                           matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
                           matrix(rnbinom(1e3, mu = 10, size = 5), ncol = 1),
                           matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1),
                           matrix(rnbinom(1e3, mu = 10, size = 5), ncol = 1),
                           matrix(rnbinom(1e3, mu = 1, size = 10), ncol = 1)))
  
  colData <- data.frame(condition = c('A','B'), replicate = c(1,1))
  rowRanges <- GenomicRanges::GRanges('chrA',
                                      IRanges::IRanges(start = seq(1,by = 500,
                                                                   length.out = nrow(countData)),width = 500))
  
  object <- epigraHMMDataSetFromMatrix(countData,colData,rowRanges = rowRanges)
  
  # Initializing
  object <- initializer(object,controlEM())
  
  # Running epigraHMM
  object <- epigraHMM(object,controlEM(maxIterEM = 5),type = 'differential',dist = 'nb')
  
  # Calling peaks with FDR thresholding
  peaks <- callPeaks(object = object,
                     hdf5 = S4Vectors::metadata(object)$output,
                     method = 0.05,
                     saveToFile = TRUE,
                     control = S4Vectors::metadata(object)$control)
  
  output <- do.call(file.path,S4Vectors::metadata(object)$control[c('tempDir','fileName')])
  output <- paste0(output,c('_peaks.bed','_prob_chrA.wig','_mixProb_EB.wig','_mixProb_BE.wig'))
  
  expect_true(all(file.exists(output)))
})
