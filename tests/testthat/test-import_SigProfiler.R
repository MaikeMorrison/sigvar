 test_that("import works", {
   expect_no_error( import_SigProfiler(system.file("extdata", "SP", package = "sigvar")))
})
