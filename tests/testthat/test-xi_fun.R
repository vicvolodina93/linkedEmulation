test_that("independent case for xi_fun works!", {
  expect_equal(xi_fun(w_i = c(1, 3, 4),
                      mu = c(0, 13, 45),
                      Sigma = diag(c(4, 10, 0.8), ncol = 3, nrow = 3),
                      params = list(delta_par = c(0.1, 1, 0.75))),
               c(0.0311863613, 0.0018656123, 0.000000000))
})

test_that("dependent case for xi_fun works!", {
  expect_equal(xi_fun(w_i = c(1, 3, 4),
                        mu = c(4, 13, 45),
                        Sigma =matrix(c(400, 0, 0, 0, 100, 39, 0, 39, 45),
                                      ncol = 3, nrow = 3),
                        params = list(delta_par = c(10, 4, 0.75)),
                        ind = "FALSE"), matrix(9.994653e-13, ncol = 1,
                                               nrow = 1))
})
