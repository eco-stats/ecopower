context("check input args")

test_that("cord object", {

  expect_error(
    equivtest(fit1.glm, effect_mat, term_cont)
  )

})

test_that("coeffs exist", {

  expect_error(
    equivtest(fit1.cord, NULL)
  )

})

test_that("both term and object0 are specified", {

  expect_error(
    equivtest(fit1.cord, effect_mat, term_cont, fit0.cord)
  )

})
