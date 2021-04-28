test_that("check epigraHMM object from matrix is created properly",{
    
    # Creating input
    mat <- matrix(rpois(4e5,10),ncol = 4)
    dsg <- data.frame(condition = c('A','A','B','B'),replicate = c(1,2,1,2))
    ranges <- GenomicRanges::GRanges('chrA',
                                     IRanges::IRanges(start = seq(1,by = 500,length.out = 1e5),
                                                      width = 500))
    
    object <- epigraHMMDataSetFromMatrix(countData = mat,colData = dsg,rowRanges = ranges)
    
    # Testing
    expect_equal(assay(object),mat,check.attributes = FALSE)
    
    expect_equal(as.matrix(assay(object,'offsets')),
                 matrix(0,nrow = 1e5,ncol = 4),check.attributes = FALSE)
    
    expect_equal(colData(object),S4Vectors::DataFrame(dsg),check.attributes = FALSE)
    
    expect_equal(rowRanges(object),ranges)
})

test_that("check epigraHMM object from matrix is sorted",{
    
    # Creating input
    mat <- matrix(rpois(4e5,10),ncol = 4)
    dsg <- data.frame(condition = c('A','B','A','B'),replicate = c(2,2,1,1))
    ranges <- GenomicRanges::GRanges('chrA',
                                     IRanges::IRanges(start = seq(1,by = 500,length.out = 1e5),
                                                      width = 500))
    
    object <- epigraHMMDataSetFromMatrix(countData = mat,colData = dsg,rowRanges = ranges)
    
    # Testing
    
    expect_equal(colData(object),S4Vectors::DataFrame(dsg[order(dsg$condition,dsg$replicate),]),
                 check.attributes = FALSE)
    
})