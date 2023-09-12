library(dplyr)
mutsig_carcinogens_relab = mutsig_carcinogens_mice_SBS %>%
  select(Sample, chemical, Tissue, mSBS1:mSBS_N3)


test_that("bootstrapping works: no groups", {
  expect_no_error(sigboot(mutsig_carcinogens_relab, K = 11, n_replicates = 10))
})

test_that("bootstrapping works: groups", {
  expect_no_error(bootstrap_fava(matrices = mutsig_carcinogens_relab,
                                 K = 11,
                                 group = c("Tissue"),
                                 seed = 1,
                                 n_replicates = 100))
  expect_warning(bootstrap_fava(matrices = mutsig_carcinogens_relab,
                                 K = 11,
                                 group = c("Tissue", "chemical"),
                                 seed = 1,
                                 n_replicates = 100))
})
