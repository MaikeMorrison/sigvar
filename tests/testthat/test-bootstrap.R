library(dplyr)

# Load data from sigvar package:
ds_info <- data(package = "sigvar")
ds_names <- ds_info$results[, "Item"]
ds_clean <- gsub(" .*$", "", ds_names)
data(list = ds_clean, package = "sigvar")


test_that("bootstrapping works for a matrix with one grouping var, multiple groups", {
  # NO WEIGHTS
  expect_no_error(sigboot(
    sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
    K = 43
  ))
  # NO K SPECIFIED WHEN EVERY COLUMN USED
  expect_identical(
    sigboot(sig_activity = ESCC_sig_activity %>% select(-c(Incidence_Level, Sample)), n_replicates = 3, group = "Country", seed = 1)$P_values,
    sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country", K = 43, seed = 1)$P_values
  )
  expect_identical(
    sigboot(
      sig_activity = ESCC_sig_activity %>% select(-c(Incidence_Level, Sample)), n_replicates = 3, group = "Country",
      seed = 1, S = ESCC_sig_similarity
    )$P_values,
    sigboot(
      sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country", K = 43,
      seed = 1, S = ESCC_sig_similarity
    )$P_values
  )
  # SIMILARITY
  expect_no_error(sigboot(
    sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
    K = 43, S = ESCC_sig_similarity
  ))
  # NORMALIZED
  expect_no_error(sigboot(
    sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
    K = 43,
    normalized = TRUE
  ))
  # NORMALIZED AND SIMILARITY
  expect_error(sigboot(
    sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
    K = 43, S = ESCC_sig_similarity,
    normalized = TRUE
  ))
})

relab_2_groups <- ESCC_sig_activity %>% dplyr::filter(Country %in% c("UK", "China"))

test_that("bootstrapping works for a matrix with one grouping var, two groups", {
  # NO WEIGHTS
  expect_no_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country", K = 43))
  # NO K SPECIFIED WHEN EVERY COLUMN USED
  expect_identical(
    sigboot(sig_activity = relab_2_groups %>% select(-c(Incidence_Level, Sample)), n_replicates = 3, group = "Country", seed = 1)$P_values,
    sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country", K = 43, seed = 1)$P_values
  )
  expect_identical(
    sigboot(
      sig_activity = relab_2_groups %>% select(-c(Incidence_Level, Sample)), n_replicates = 3, group = "Country",
      seed = 1, S = ESCC_sig_similarity
    )$P_values,
    sigboot(
      sig_activity = relab_2_groups, n_replicates = 3, group = "Country", K = 43,
      seed = 1, S = ESCC_sig_similarity
    )$P_values
  )
  # SIMILARITY
  expect_no_error(sigboot(
    sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
    K = 43, S = ESCC_sig_similarity
  ))
  # NORMALIZED
  expect_no_error(sigboot(
    sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
    K = 43,
    normalized = TRUE
  ))
  # NORMALIZED AND SIMILARITY
  expect_error(sigboot(
    sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
    K = 43, S = ESCC_sig_similarity,
    normalized = TRUE
  ))
})


test_that("bootstrapping yields expected error for one matrix", {
  expect_error(sigboot(
    sig_activity = dplyr::filter(ESCC_sig_activity, Country == "China"),
    n_replicates = 3, K = 43, S = ESCC_sig_similarity, group = "Country"
  ))
})
