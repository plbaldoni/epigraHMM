test_that("check if addOffsets gives correct output",{
     object <- epigraHMMDataSetFromMatrix(matrix(rpois(4e5,10),ncol = 4),
                                          data.frame(condition = c('A','A','B','B'), 
                                                     replicate = c(1,2,1,2)))
     
     # Checking input
     expect_error(addOffsets(object,matrix(NA,ncol = 4,nrow = 1e5)))
     expect_error(addOffsets(object,matrix(Inf,ncol = 4,nrow = 1e5)))
     expect_error(addOffsets(object,matrix(NaN,ncol = 4,nrow = 1e5)))
     
     # Checking output
     obj <- assay(addOffsets(addOffsets(object,matrix(pi,ncol = 4,nrow = 1e5)),
                             matrix(pi,ncol = 4,nrow = 1e5)),'offsets')
     expect_true(all(as.matrix(obj)==matrix(2*pi,nrow = 1e5,ncol = 4)))
})