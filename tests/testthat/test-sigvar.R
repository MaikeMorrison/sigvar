S = cossim(ref_sigs = as.matrix(Sherlock_LCINS_SBS.refs[,1:14]))

sigvar_LCINS = sigvar(sig_activity = zhang_sig_activity, K = 14,
       group = "passive_smoking", S = S)

test_that("sigvar works", {
  expect_equal(sigvar_LCINS[1,2] , 0.1178677128080407310318)
  expect_equal(sigvar_LCINS[1,3] , 0.3260461117384914619954)
  expect_equal(sigvar_LCINS[2,2] , 0.130207605808201842823)
  expect_equal(sigvar_LCINS[2,3] , 0.3496011620113261830767)
})
