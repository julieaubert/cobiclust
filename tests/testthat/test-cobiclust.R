test_that("arguments for cobiclust main function", {

  ## parameters
  expect_is(x, "matrix")
  expect_true(K >= 1)
  expect_true(K < nrow(x))
  expect_true(G >= 1)
  expect_true(K < ncol(x))
})
