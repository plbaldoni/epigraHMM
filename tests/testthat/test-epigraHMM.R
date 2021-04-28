test_that("check if HDF5 output has the correct attributes (consensus peak calling with NB distribution)",{
    
    # Creating dummy object
    set.seed(210423)
    countData <- rbind(matrix(rnbinom(2e3,mu = 2,size = 10),ncol = 2),
                       matrix(rnbinom(4e3,mu = 10,size = 2),ncol = 2),
                       matrix(rnbinom(2e3,mu = 2,size = 10),ncol = 2))
    
    colData <- data.frame(condition = c('A','A'), replicate = c(1,2))
    
    rowRanges <- GenomicRanges::GRanges('chrA',
                                        IRanges::IRanges(start = seq(from = 1, length.out = 4e3,by = 250),width = 250))
    
    object <- epigraHMMDataSetFromMatrix(countData,colData,rowRanges)
    
    # Initializing
    object <- initializer(object,controlEM())
    
    # Running epigraHMM
    object <- epigraHMM(object,controlEM(),type = 'consensus',dist = 'nb')
    
    # Tests
    output <- metadata(object)$output
    
    ## File exists
    expect_true(file.exists(output))
    
    ## Check names
    expect_true(all.equal(rhdf5::h5ls(output)$name,c("logBP","logFP","logLikelihood",
                                                     "logProb1","logProb2","viterbi")))
    
    ## Check dimensions
    expect_true(all(unlist(lapply(rhdf5::h5ls(output)$name,function(x){
        nrow(rhdf5::h5read(output,x)) == 4e3
    }))))
    
    ## Check content
    expect_true(all(unlist(lapply(rhdf5::h5ls(output)$name,function(x){
        content <- rhdf5::h5read(output,x)
        sum(is.na(content) | is.nan(content)) == 0
    }))))
})

test_that("check if HDF5 output has the correct attributes (consensus peak calling with ZINB distribution)",{
    
    # Creating dummy object
    set.seed(210423)
    countData <- rbind(matrix(rnbinom(2e3,mu = 2,size = 10),ncol = 2),
                       matrix(rnbinom(4e3,mu = 10,size = 2),ncol = 2),
                       matrix(rnbinom(2e3,mu = 2,size = 10),ncol = 2))
    
    colData <- data.frame(condition = c('A','A'), replicate = c(1,2))
    
    rowRanges <- GenomicRanges::GRanges('chrA',
                                        IRanges::IRanges(start = seq(from = 1, length.out = 4e3,by = 250),width = 250))
    
    object <- epigraHMMDataSetFromMatrix(countData,colData,rowRanges)
    
    # Initializing
    object <- initializer(object,controlEM())
    
    # Running epigraHMM
    object <- epigraHMM(object,controlEM(),type = 'consensus',dist = 'zinb')
    
    # Tests
    output <- metadata(object)$output
    
    ## File exists
    expect_true(file.exists(output))
    
    ## Check names
    expect_true(all.equal(rhdf5::h5ls(output)$name,c("logBP","logFP","logLikelihood",
                                                     "logProb1","logProb2","viterbi")))
    
    ## Check dimensions
    expect_true(all(unlist(lapply(rhdf5::h5ls(output)$name,function(x){
        nrow(rhdf5::h5read(output,x)) == 4e3
    }))))
    
    ## Check content
    expect_true(all(unlist(lapply(rhdf5::h5ls(output)$name,function(x){
        content <- rhdf5::h5read(output,x)
        sum(is.na(content) | is.nan(content)) == 0
    }))))
})

test_that("check if HDF5 output has the correct attributes (differential peak with NB distribution)",{
    
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
    object <- epigraHMM(object,controlEM(),type = 'differential',dist = 'nb')
    
    # Tests
    output <- metadata(object)$output
    
    ## File exists
    expect_true(file.exists(output))
    
    ## Check names
    expect_true(all.equal(rhdf5::h5ls(output)$name,c("logBP","logFP","logLikelihood",
                                                     "logProb1","logProb2","mixturePatterns",
                                                     "mixtureProb","viterbi")))
    
    ## Check dimensions
    expect_true(all(unlist(lapply(rhdf5::h5ls(output)$name[-6],function(x){
        nrow(rhdf5::h5read(output,x)) == 7e3
    }))))
    
    ## Check content
    expect_true(all(unlist(lapply(rhdf5::h5ls(output)$name,function(x){
        content <- rhdf5::h5read(output,x)
        sum(is.na(content) | is.nan(content)) == 0
    }))))
})

test_that("check if HDF5 output has the correct attributes (differential peak with ZINB distribution)",{
    
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
    object <- epigraHMM(object,controlEM(),type = 'differential',dist = 'zinb')
    
    # Tests
    output <- metadata(object)$output
    
    ## File exists
    expect_true(file.exists(output))
    
    ## Check names
    expect_true(all.equal(rhdf5::h5ls(output)$name,c("logBP","logFP","logLikelihood",
                                                     "logProb1","logProb2","mixturePatterns",
                                                     "mixtureProb","viterbi")))
    
    ## Check dimensions
    expect_true(all(unlist(lapply(rhdf5::h5ls(output)$name[-6],function(x){
        nrow(rhdf5::h5read(output,x)) == 7e3
    }))))
    
    ## Check content
    expect_true(all(unlist(lapply(rhdf5::h5ls(output)$name,function(x){
        content <- rhdf5::h5read(output,x)
        sum(is.na(content) | is.nan(content)) == 0
    }))))
})