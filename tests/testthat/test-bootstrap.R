library(dplyr)

# sig_activity = ESCC_sig_activity %>% dplyr::filter(Country %in% c("UK", "China", "Japan"))
# n_replicates=100
# K = 43
# group = "Country"# c("Country", "Incidence_Level")
# S = NULL
# normalized = FALSE
# seed = NULL
# save_replicates = FALSE
# alternative = "less"
# w = NULL
# time = NULL



test_that("bootstrapping works for a matrix with one grouping var, multiple groups", {
  # NO WEIGHTS
  expect_no_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
                                K = 43))
  # NO K SPECIFIED WHEN EVERY COLUMN USED
  expect_identical(sigboot(sig_activity = ESCC_sig_activity %>% select(-c(Incidence_Level, Sample)), n_replicates = 3, group = "Country", seed = 1)$P_values,
                   sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",K = 43, seed = 1)$P_values)
  expect_identical(sigboot(sig_activity = ESCC_sig_activity %>% select(-c(Incidence_Level, Sample)), n_replicates = 3, group = "Country",
                                  seed = 1, S = ESCC_sig_similarity)$P_values,
                   sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",K = 43,
                                  seed = 1, S = ESCC_sig_similarity)$P_values)
  # SIMILARITY
  expect_no_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
                                K = 43, S = ESCC_sig_similarity))
  # # TIME
  # expect_no_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
  #                               K = 43))
  # # W
  # expect_no_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
                                # K = 43, w = time_weights(times = sig_activity$timepoint, group = sig_activity$subject)))
  # CONFIRM TIME AND W IDENTICAL w same seed
  # expect_identical(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 4, group = "Country",
  #                                K = 43, seed = 1)$P_values,
  #                  sigboot(sig_activity = ESCC_sig_activity, n_replicates = 4, group = "Country",
  #                                K = 43, w = time_weights(times = sig_activity$timepoint, group = sig_activity$subject), seed = 1)$P_values)
  # # SIMILARITY AND TIME
  # expect_no_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
  #                               K = 43, S = ESCC_sig_similarity))
  # # SIMILARITY AND W
  # expect_no_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
  #                               K = 43, w = time_weights(times = sig_activity$timepoint, group = sig_activity$subject),
  #                                S = ESCC_sig_similarity))
  # # TIME AND W
  # expect_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
  #                            K = 43, w = time_weights(times = sig_activity$timepoint, group = sig_activity$subject),
  #                             time = "timepoint"))

  # NORMALIZED
  expect_no_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
                                K = 43,
                                 normalized = TRUE))
  # NORMALIZED AND SIMILARITY
  expect_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
                             K = 43, S = ESCC_sig_similarity,
                              normalized = TRUE))
  # # NORMALIZED AND TIME
  # expect_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
  #                            K = 43,time = "timepoint",
  #                             normalized = TRUE))
  # # NORMALIZED AND W
  # expect_error(sigboot(sig_activity = ESCC_sig_activity, n_replicates = 3, group = "Country",
  #                            K = 43, w = time_weights(times = sig_activity$timepoint, group = sig_activity$subject),
  #                             normalized = TRUE))

})

relab_2_groups = ESCC_sig_activity %>% dplyr::filter(Country %in% c("UK", "China"))

test_that("bootstrapping works for a matrix with one grouping var, two groups", {
  # NO WEIGHTS
  expect_no_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",K = 43))
  # NO K SPECIFIED WHEN EVERY COLUMN USED
  expect_identical(sigboot(sig_activity = relab_2_groups %>% select(-c(Incidence_Level, Sample)), n_replicates = 3, group = "Country", seed = 1)$P_values,
                   sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",K = 43, seed = 1)$P_values)
  expect_identical(sigboot(sig_activity = relab_2_groups %>% select(-c(Incidence_Level, Sample)), n_replicates = 3, group = "Country",
                                  seed = 1, S = ESCC_sig_similarity)$P_values,
                   sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",K = 43,
                                  seed = 1, S = ESCC_sig_similarity)$P_values)
  # SIMILARITY
  expect_no_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
                                K = 43, S = ESCC_sig_similarity))
  # # TIME
  # expect_no_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
  #                               K = 43))
  # # W
  # expect_no_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
  #                               K = 43, w = time_weights(times = relab_2_groups$timepoint, group = relab_2_groups$subject)))
  # # CONFIRM TIME AND W IDENTICAL w same seed
  # expect_identical(sigboot(sig_activity = relab_2_groups, n_replicates = 4, group = "Country",
  #                                K = 43, seed = 1)$P_values,
  #                  sigboot(sig_activity = relab_2_groups, n_replicates = 4, group = "Country",
  #                                K = 43, w = time_weights(times = relab_2_groups$timepoint, group = relab_2_groups$subject), seed = 1)$P_values)
  # # SIMILARITY AND TIME
  # expect_no_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
  #                               K = 43, S = ESCC_sig_similarity))
  # # SIMILARITY AND W
  # expect_no_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
  #                               K = 43, w = time_weights(times = relab_2_groups$timepoint, group = relab_2_groups$subject),
  #                                S = ESCC_sig_similarity))
  # # TIME AND W
  # expect_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
  #                            K = 43, w = time_weights(times = relab_2_groups$timepoint, group = relab_2_groups$subject),
  #                             time = "timepoint"))
  # NORMALIZED
  expect_no_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
                                K = 43,
                                 normalized = TRUE))
  # NORMALIZED AND SIMILARITY
  expect_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
                             K = 43, S = ESCC_sig_similarity,
                              normalized = TRUE))
  # # NORMALIZED AND TIME
  # expect_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
  #                            K = 43,time = "timepoint",
  #                             normalized = TRUE))
  # # NORMALIZED AND W
  # expect_error(sigboot(sig_activity = relab_2_groups, n_replicates = 3, group = "Country",
  #                            K = 43, w = time_weights(times = relab_2_groups$timepoint, group = relab_2_groups$subject),
  #                             normalized = TRUE))

})



test_that("bootstrapping yields expected error for one matrix", {
  expect_error(sigboot(sig_activity = dplyr::filter(ESCC_sig_activity, Country == "China"),
                              n_replicates = 3,K = 43, S = ESCC_sig_similarity, group = "Country"))
})


# library(dplyr)
# test_groups = sig_activity %>%
#   mutate(Abx = ifelse(timepoint < 29, "Before", ifelse(timepoint > 34, "After", "During")),
#          .before = 1) %>% filter(Abx != "During")
#
#
# test_that("bootstrapping works for a matrix with two grouping vars, multiple groups", {
#   # NO WEIGHTS
#   expect_no_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                                 K = 43))
#   # NO K SPECIFIED WHEN EVERY COLUMN USED
#   expect_identical(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"), seed = 1)$P_values,
#                    sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),K = 43, seed = 1)$P_values)
#   expect_identical(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                                   seed = 1, S = ESCC_sig_similarity)$P_values,
#                    sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),K = 43,
#                                   seed = 1, S = ESCC_sig_similarity)$P_values)
#   # SIMILARITY
#   expect_no_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                                 K = 43, S = ESCC_sig_similarity))
#   # TIME
#   expect_no_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                                 K = 43))
#   # W
#   expect_no_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                                 K = 43, w = time_weights(times = test_groups$timepoint, group = test_groups$subject)))
#   # CONFIRM TIME AND W IDENTICAL w same seed
#   expect_identical(sigboot(sig_activity = test_groups, n_replicates = 4, group = c("subject", "Abx"),
#                                  K = 43, seed = 1)$P_values,
#                    sigboot(sig_activity = test_groups, n_replicates = 4, group = c("subject", "Abx"),
#                                  K = 43, w = time_weights(times = test_groups$timepoint, group = paste0(test_groups$Abx, test_groups$subject)), seed = 1)$P_values)
#   # SIMILARITY AND TIME
#   expect_no_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                                 K = 43, S = ESCC_sig_similarity))
#   # SIMILARITY AND W
#   expect_no_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                                 K = 43, w = time_weights(times = test_groups$timepoint, group = test_groups$subject),
#                                  S = ESCC_sig_similarity))
#   # TIME AND W
#   expect_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                              K = 43, w = time_weights(times = test_groups$timepoint, group = test_groups$subject),
#                               time = "timepoint"))
#   # NORMALIZED
#   expect_no_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                                 K = 43,
#                                  normalized = TRUE))
#   # NORMALIZED AND SIMILARITY
#   expect_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                              K = 43, S = ESCC_sig_similarity,
#                               normalized = TRUE))
#   # NORMALIZED AND TIME
#   expect_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                              K = 43,time = "timepoint",
#                               normalized = TRUE))
#   # NORMALIZED AND W
#   expect_error(sigboot(sig_activity = test_groups, n_replicates = 3, group = c("subject", "Abx"),
#                              K = 43, w = time_weights(times = test_groups$timepoint, group = test_groups$subject),
#                               normalized = TRUE))
#
# })

