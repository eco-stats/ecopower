context("check input args")

test_that("cord object", {

  expect_error(
    powersim(fit1.glm, effect_mat, term_cont)
  )

})

test_that("coeffs exist", {

  expect_error(
    powersim(fit1.cord, NULL, term_cont)
  )

})
