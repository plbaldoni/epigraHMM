test_that("check dot/dash in path",{
        
        tmp <- file.path(tempdir(),'dot.dash-path')
        
        file <- tempfile(pattern = 'epigraHMM',tmpdir = tmp,fileext = '.h5')
        expect_true(checkPath(file) == file)
})

test_that("check dot/dash in path with existing file",{
        
        tmp <- file.path(tempdir(),'dot.dash-path')
        
        file <- tempfile(pattern = 'epigraHMM',tmpdir = tmp,fileext = '.h5')
        
        dir.create(tmp)
        dump <- file.create(file,showWarnings = FALSE)
        
        # Should return a path with a new file
        expect_false(file.exists(checkPath(file)))
})
