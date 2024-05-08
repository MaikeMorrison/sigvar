SBS_table_test = sigvar::COSMIC3.3.1_SBS[,1:5]
SBS_table_big = SBS_table_test[,-1]*20

sigvar::plot_SBS_spectrum(sigvar::COSMIC3.3.1_SBS[,2:4])


test_that("plot_SBS_spectrum checks for non numeric cols", {
  expect_warning(plot_SBS_spectrum(SBS_table_test))
})


test_that("plot_SBS_spectrum checks for 96 rows", {
  expect_error(plot_SBS_spectrum(SBS_table_test[1:10,]))
})


test_that("plot_SBS_spectrum checks for colsums to 1", {
  expect_warning(plot_SBS_spectrum(SBS_table_big))
})


test_that("plot_SBS_spectrum works fine when data is correctly formatted", {
  expect_no_warning(plot_SBS_spectrum(SBS_table_test[,-1]))
})
#
