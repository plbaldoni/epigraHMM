test_that("check normalizeCounts w/ dummy data",{
    
    # Creating input
    mat <- matrix(rep(100,4e5),ncol = 4)
    dsg <- data.frame(condition = c('A','A','B','B'),replicate = c(1,2,1,2))
    ranges <- GenomicRanges::GRanges('chrA',
                                     IRanges::IRanges(start = seq(1,by = 500,length.out = 1e5),
                                                      width = 500))
    
    object <- epigraHMMDataSetFromMatrix(countData = mat,colData = dsg,rowRanges = ranges)
    
    objectNormalized <- normalizeCounts(object,controlEM())
    
    expect_equal(unique(c(as.matrix(assay(objectNormalized,'offsets')))),0)
    
    # Creating input
    mat <- matrix(c(1,2,1,2),ncol = 4,nrow = 1e5,byrow = TRUE)

    object <- epigraHMMDataSetFromMatrix(countData = mat,colData = dsg,rowRanges = ranges)
    
    objectNormalized <- normalizeCounts(object,controlEM())
    
    expect_true(all(Matrix::rowSums(assay(objectNormalized,'offsets')) == 0))
})