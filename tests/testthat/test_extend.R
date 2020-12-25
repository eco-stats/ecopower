context("check input args")

test_that("cord object", {

  expect_error(
    extend(fit1.glm)
  )

})

context("compute_N")

test_that("n_replicate overwrites N", {

  expect_warning(
    extend(fit1.cord, n_replicate=2)
  )

})

context("check output data.frame")

test_that("N without newdata", {

  for (fit in list(fit1.cord, fit2.cord, fit_fac4.cord, fit_mth.cord, fit_mix.cord)) {

    returned = extend(fit, N=nrow(spiddat)*2)$data
    expect_equal(nrow(returned), 56)
    expect_equal(ncol(returned), ncol(spiddat)+ncol(X))
    expect_equal(
      colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
    )
    expect_equal(
      colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(X))],
      colnames(X)
    )
  }

})

test_that("n_replicate without newdata", {

  for (fit in list(fit1.cord, fit2.cord, fit_fac4.cord, fit_mth.cord, fit_mix.cord)) {

    returned = extend(fit, N=NULL, n_replicate=3)$data
    expect_equal(nrow(returned), 84)
    expect_equal(ncol(returned), ncol(spiddat)+ncol(X))
    expect_equal(
      colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
    )
    expect_equal(
      colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(X))],
      colnames(X)
    )
  }

})

test_that("N with newdata", {

  for (fit in list(fit_fac2.cord, fit_mth.cord)) {

    returned = extend(fit, N=30, newdata=X_new)$data
    expect_equal(nrow(returned), 30)
    expect_equal(ncol(returned), ncol(spiddat)+ncol(X))
    expect_equal(
      colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
    )
    expect_equal(
      colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(X))],
      colnames(X)
    )
  }

})

test_that("n_replicate with newdata", {

  for (fit in list(fit_fac2.cord, fit_mth.cord)) {

    returned = extend(fit, N=NULL, newdata=X_new, n_replicate=2)$data
    expect_equal(nrow(returned), 56)
    expect_equal(ncol(returned), ncol(spiddat)+ncol(X))
    expect_equal(
      colnames(returned)[1:ncol(spiddat)], colnames(spiddat)
    )
    expect_equal(
      colnames(returned)[(ncol(spiddat)+1):(ncol(spiddat)+ncol(X))],
      colnames(X)
    )
  }

})

test_that("intercept model", {

  returned = extend(fit0.cord, N=nrow(spiddat)*2)$data
  returned1 = extend(fit0.cord, N=NULL, n_replicate=2)$data
  expect_equal(nrow(returned), 56)
  expect_equal(nrow(returned), nrow(returned1))
  expect_equal(colnames(returned), colnames(returned1))

})

context("get_new_fit")

test_that("do.fit=TRUE returns manyglm object", {

  returned = extend(fit1.cord, do.fit=TRUE)
  expect_is(returned, "manyglm")

})

test_that("do.fit=FALSE returns a data.frame", {

  returned = extend(fit1.cord)$data
  expect_is(returned, "data.frame")

})
