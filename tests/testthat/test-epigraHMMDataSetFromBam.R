test_that("check output from epigraHMM object (bam)",{
    
    if ("chromstaRData" %in% rownames(installed.packages())) {
        bamFiles <- system.file("extdata","euratrans",
                                "lv-H3K27me3-SHR-male-bio2-tech1.bam",
                                package="chromstaRData")   
        
        colData <- data.frame(condition = 'SHR', replicate = 1)
        
        object <- epigraHMMDataSetFromBam(bamFiles = bamFiles,
                                          colData = colData,
                                          genome = 'rn4',
                                          windowSize = 10000,
                                          gapTrack = TRUE,
                                          blackList = TRUE)
        # Check genome
        expect_true(all(genome(object)=='rn4'))
        
        # Check window size
        expect_true(sum(width(object)==10000)>=(nrow(object)-1))
        
        # Check genome via GRanges
        objectOne <- epigraHMMDataSetFromBam(bamFiles = bamFiles,
                                             colData = colData,
                                             genome = rowRanges(object)[1],
                                             windowSize = 10000,
                                             gapTrack = TRUE,
                                             blackList = TRUE)
        
        expect_equal(assay(objectOne)[1],assay(object)[1])
        expect_equal(width(objectOne),10000)
        
        # Check discards
        objectTwo <- epigraHMMDataSetFromBam(bamFiles = bamFiles,
                                             colData = colData,
                                             genome = rowRanges(object)[1:3],
                                             windowSize = 10000,
                                             gapTrack = rowRanges(object)[2],
                                             blackList = rowRanges(object)[3])
        
        expect_equal(assay(objectTwo)[1],assay(objectOne)[1])
        expect_equal(width(objectTwo),10000)
    }
})