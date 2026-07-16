SPfolder <- system.file("extdata", "SP", package = "sigvar")

A <- list(
  "DBS78_De-Novo_Solution"=
  tibble("Samples"=paste0("Sample",seq_len(10)),
         "DBS78A"=c(12,  0,  9, 22, 20, 19, 15, 27,  5, 21),
         "DBS78B"=c(6,  0,  1,  6,  3,  0, 10,  5,  6,  6)
  ))

A_error = A
A_error[[1]][2,3] = 7

test_that("import_SigProfiler works", {
  expect_no_warning(import_SigProfiler(SPfolder))
  expect_true(all(import_SigProfiler(SPfolder)$`DBS78_De-Novo_Solution`==A$`DBS78_De-Novo_Solution`) )
  expect_error(import_SigProfiler("test"))
  expect_false( all(import_SigProfiler(SPfolder)$`DBS78_De-Novo_Solution`==A_error$`DBS78_De-Novo_Solution`) )
})
