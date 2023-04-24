test_that("import works", {
  expect_equal( import_SigProfiler(system.file("extdata", "example_SigProfiler_results", package = "sigFAVA")), internal_for_import_test)
})
