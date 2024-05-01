sig_activity = ESCC_sig_activity
group = "Country"
K=ncol(sig_activity)-3
max_dotsize = 5
pivot = FALSE
median = FALSE
normalized = FALSE
threshold = 0
facet = "Incidence_Level"


# WITH FACET
test_that("mean normalized works with facet", {
  expect_warning(plot_dots(sig_activity = sig_activity, group = "Country",
                           facet = "Incidence_Level", K = 43))
  expect_no_warning(plot_dots(sig_activity = ESCC_sig_activity, group = "Country",
                              facet = "Incidence_Level", K = 43))
})


test_that("mean unnormalized works with facet", {
  expect_no_error(plot_dots(sig_activity = sig_activity, group = "Country",
                            facet = "Incidence_Level", normalized = FALSE, K = 43))
  expect_no_warning(plot_dots(sig_activity = sig_activity, group = "Country",
                              facet = "Incidence_Level", normalized = FALSE, K = 43))
})

test_that("median normalized works with facet", {
  expect_warning(plot_dots(sig_activity = sig_activity, median = TRUE, group = "Country",
                           facet = "Incidence_Level", K = 43))
  expect_no_warning(plot_dots(sig_activity = ESCC_sig_activity, median = TRUE, group = "Country",
                              facet = "Incidence_Level", K = 43))
})


test_that("median unnormalized works with facet", {
  expect_no_error(plot_dots(sig_activity = sig_activity, median = TRUE, group = "Country",
                            facet = "Incidence_Level", normalized = FALSE, K = 43))
  expect_no_warning(plot_dots(sig_activity = sig_activity, median = TRUE, group = "Country",
                              facet = "Incidence_Level", normalized = FALSE, K = 43))
})



# NO FACET
test_that("mean normalized works without facet", {
  expect_warning(plot_dots(sig_activity = sig_activity, group = "Country",
                           K = 43))
  expect_no_warning(plot_dots(sig_activity = ESCC_sig_activity, group = "Country",
                              K = 43))
})


test_that("mean unnormalized works without facet", {
  expect_no_error(plot_dots(sig_activity = sig_activity, group = "Country",
                            normalized = FALSE, K = 43))
  expect_no_warning(plot_dots(sig_activity = sig_activity, group = "Country",
                              normalized = FALSE, K = 43))
})

test_that("median normalized works without facet", {
  expect_warning(plot_dots(sig_activity = sig_activity, median = TRUE, group = "Country",
                           K = 43))
  expect_no_warning(plot_dots(sig_activity = ESCC_sig_activity, median = TRUE, group = "Country",
                              K = 43))
})


test_that("median unnormalized works without facet", {
  expect_no_error(plot_dots(sig_activity = sig_activity, median = TRUE, group = "Country",
                            normalized = FALSE, K = 43))
  expect_no_warning(plot_dots(sig_activity = sig_activity, median = TRUE, group = "Country",
                              normalized = FALSE, K = 43))
})

