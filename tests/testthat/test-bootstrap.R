library(dplyr)

sig_activity = ESCC_sig_activity %>% filter(Country %in% c("UK", "China"))
n_replicates=100
K = 43
group = "Country"
S = NULL
normalized = FALSE
seed = NULL
save_replicates = FALSE
alternative = "lesser"

test_that("bootstrapping works with just two groups", {
  expect_no_error(sig_boot(sig_activity = smoker_sigs_chen, n_replicates = 10, K = 3, group = "Smoker"))
})




# OLD
# mutsig_carcinogens_relab = mutsig_carcinogens_mice_SBS %>%
# select(Sample, chemical, Tissue, mSBS1:mSBS_N3)
# test_that("bootstrapping works: no groups", {
#   expect_no_error(sigboot(mutsig_carcinogens_relab, K = 11, n_replicates = 10))
#   expect_no_error(sigboot(ESCC_sig_activity, K = 43, n_replicates = 10))
# })
#
# test_that("bootstrapping works: groups", {
#   expect_no_error(bootstrap_fava(matrices = mutsig_carcinogens_relab,
#                                  K = 11,
#                                  group = c("Tissue"),
#                                  seed = 1,
#                                  n_replicates = 100))
#
#   expect_no_error(sigboot(ESCC_sig_activity, group = "Country", K = 43, n_replicates = 10))
#
#   expect_warning(bootstrap_fava(matrices = mutsig_carcinogens_relab,
#                                  K = 11,
#                                  group = c("Tissue", "chemical"),
#                                  seed = 1,
#                                  n_replicates = 100))
# })
#
# test_that("bootstrapping works: similarity and group", {
#   expect_no_error(sigboot(ESCC_sig_activity, group = "Country", K = 43, n_replicates = 10, S = ESCC_sig_similarity))
# })
