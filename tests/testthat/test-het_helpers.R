q1 <- c(1, 0, 0, 0)
q2 <- c(0.5, 0.5, 0, 0)
q3 <- c(1 / 4, 1 / 4, 1 / 4, 1 / 4)
q4 <- c(0, 0, 1, 0)

relative_abundances <- matrix(c(q1, q2, q3, q4),
  byrow = TRUE, nrow = 4
)

test_that("het works", {
  expect_equal(het(c(0.5, 0.5)), 0.5)
  expect_equal(mean(vapply(X = list(q1, q2, q3, q4), FUN = het, FUN.VALUE = numeric(1))), 0.3125)
})


test_that("het_mean works", {
  expect_equal(het_mean(relative_abundances), 0.3125)
})

# Compute mean Gini-Simpson index ignoring
# rows 2 and 3
row_weights <- c(0.5, 0, 0, 0.5)
het_mean(relative_abundances, w = row_weights)

# Compute mean Gini-Simpson index assuming that
# categories 1 and 2 are identical:
similarity_matrix <- diag(4)
similarity_matrix[1, 2] <- 1
similarity_matrix[2, 1] <- 1

# Assume categories 1 and 2 are identical AND
# ignore rows 2 and 4:
row_weights <- c(0.5, 0, 0.5, 0)


test_that("weighted het_mean works", {
  expect_equal(het_mean(relative_abundances, S = similarity_matrix), 0.15625)
  expect_equal(het_mean(relative_abundances, w = row_weights, S = similarity_matrix), 0.3125)
})

test_that("het_pooled works", {
  expect_equal(het_pooled(relative_abundances), 0.671875)
})

test_that("weighted het_pooled works", {
  expect_equal(het_pooled(relative_abundances, S = similarity_matrix), 0.5078125)
  expect_equal(het_pooled(relative_abundances, w = row_weights, S = similarity_matrix), 0.40625)
})


test_that("fst_norm works", {
  expect_no_error(fst_norm(relative_abundances))
})


test_that("time_weights works", {
  expect_no_error(time_weights(times = c(
    1, 8, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31,
    32, 33, 34, 35, 36, 37, 38, 39, 44, 50, 57, 64
  )))
})

test_that("fst works", {
  expect_no_error(fst(relative_abundances))
  expect_equal(fst(relative_abundances, w = row_weights), 1 / 3)
  expect_no_error(fst(relative_abundances, w = row_weights, S = similarity_matrix))
})
