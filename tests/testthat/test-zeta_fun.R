test_that("independent case works", {
  expect_equal(zeta_fun(w_i = c(1, 3, 4), w_j = c(5, 0.5, 0.9),
                        mu = c(0, 13, 45), Sigma = diag(c(4, 10, 0.8),
                                                        ncol = 3, nrow = 3),
                        params = list(delta_par = c(0.1, 1, 0.75))),
               c(0.000000e+00, 0.0000142954891, 0.000000e+00))
})

test_that("dependent case works", {
  expect_equal(zeta_fun(w_i = c(1, 3, 4), w_j = c(5, 0.5, 0.9),
                        mu = c(4, 13, 45),
                        Sigma =matrix(c(400, 0, 0, 0, 100, 39, 0, 39, 45),
                                      ncol = 3, nrow = 3),
                        params = list(delta_par = c(10, 4, 0.75)),
                        ind = "FALSE"), matrix(9.110628e-18, ncol = 1,
                                               nrow = 1))
})
