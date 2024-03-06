
test_that("equivtest - class of object: cord object", {

  expect_error(
    equivtest(fit1.glm, effect_mat, term_cont)
  )

})


test_that("equivtest - check input args: coeffs exist", {

  expect_error(
    equivtest(fit1.cord, NULL)
  )

})

test_that("both term and object0 are specified", {

  expect_error(
    equivtest(fit1.cord, effect_mat, term_cont, fit0.cord)
  )

})
