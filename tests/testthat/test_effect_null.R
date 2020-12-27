context("effect_null - class of object")

test_that("manyglm object", {

  expect_error(
    effect_null(fit1.cord, term_cont)
  )

})

context("effect_null matrix")

test_that("term_cont", {

  for (fit in list(fit1.glm, fit2.glm)) {
    returned = effect_null(fit, term_cont)

    expect_equal(nrow(returned), nrow(fit$coefficients))
    expect_equal(ncol(returned), ncol(fit$coefficients))
    expect_equal(rownames(returned), rownames(fit$coefficients))
    expect_equal(colnames(returned), colnames(fit$coefficients))
    expect_equal(sum(returned[term_cont,]), 0)
  }

})

test_that("term_factors", {

  for (fit in list(fit_fac2.glm, fit_fac4.glm)) {

    term_fac = labels(terms(fit))
    returned = effect_null(fit, term_fac)

    expect_equal(nrow(returned), nrow(fit$coefficients))
    expect_equal(ncol(returned), ncol(fit$coefficients))
    expect_equal(colnames(returned), colnames(fit$coefficients))
    expect_equal(sum(returned[2:nrow(returned),]), 0)

    for (rowName in rownames(returned)[2:nrow(returned)]) {
      expect_true(startsWith(rowName, term_fac))
    }
  }

})
