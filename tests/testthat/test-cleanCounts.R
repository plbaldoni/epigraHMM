test_that("check output from cleanCounts",{
    
     # Creating dummy object
     gc <- rbeta(3e3,50,50)
     
     countData <- list('counts' = rbind(matrix(rnbinom(2e3,mu = 7.5,size = 10),ncol = 1),
                                        matrix(rnbinom(3e3,mu = exp(0.5 + 8*gc),size = 5),ncol = 1),
                                        matrix(rnbinom(2e3,mu = 7.5,size = 10),ncol = 1)),
                       'gc' = matrix(c(rbeta(2e3,50,50),gc,rbeta(2e3,50,50)),ncol = 1))
     
     colData <- data.frame(condition = 'A', replicate = 1)
     object <- epigraHMMDataSetFromMatrix(countData,colData)
     
     # Initializing
     object <- initializer(object = object,controlEM())
     
     # Cleaning counts
     object <- cleanCounts(object = object,effectNames = 'gc',byNames = 'peaks')
     
     # Check dimensions
     
     output <- as.matrix(assay(object,'offsets'))
         
     expect_equal(dim(output),dim(countData$counts))
     
     # Check values
     expect_false(any(is.na(output) | is.nan(output) | is.infinite(output)))
})
