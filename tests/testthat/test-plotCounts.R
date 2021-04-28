test_that("check if plot is working as it should",{
    
    # Creating input
    mat <- matrix(rpois(4e5,10),ncol = 4)
    dsg <- data.frame(condition = c('A','A','B','B'),replicate = c(1,2,1,2))
    ranges <- GenomicRanges::GRanges('chrA',
                                     IRanges::IRanges(start = seq(1,by = 500,length.out = 1e5),
                                                      width = 500))
    
    object <- epigraHMMDataSetFromMatrix(countData = mat,colData = dsg,rowRanges = ranges)
    
    fig <- plotCounts(object,c(100,200))
    expect_equal(methods::is(fig),'gg')
})