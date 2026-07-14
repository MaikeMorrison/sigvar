A <- matrix(
  c(
    .5, .3, .2,
    .4, .2, .4,
    .5, .4, .1,
    .6, .1, .3,
    .2, 0, .8
  ),
  ncol = 5
)
A_error = A
A_error[2,3] = 7

test_that("cossim works", {
  expect_no_warning(cossim(A))
  expect_all_equal(diag(cossim(A)), 1)
  expect_error(cossim("test"))
  expect_error(cossim(A_error))
})
