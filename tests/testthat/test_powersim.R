test_that("powersim - class of object: cord object", {

  expect_error(
    powersim(fit1.glm, effect_mat, term_cont)
  )

})

test_that("powersim - check input args: coeffs exist", {

  expect_error(
    powersim(fit1.cord, NULL, term_cont)
  )

})
