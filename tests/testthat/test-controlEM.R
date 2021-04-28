test_that("check temporary directory",{
    expect_error(controlEM(tempDir = file.path(tempdir(),'dirdoesntexist')))
    expect_true(dir.exists(controlEM()$tempDir))
})
