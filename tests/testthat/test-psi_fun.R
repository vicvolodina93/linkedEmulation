test_that("independent case for psi_fun works!", {
  expect_equal(psi_fun(w_j = c(1, 3, 4),
                       mu = c(0, 13, 45),
                       Sigma = diag(c(4, 10, 0.8), ncol = 3, nrow = 3),
                       params = list(delta_par = c(0.1, 1, 0.75))),
               c(0.044021643, 0.008960624, 0.000000000))
})

test_that("dependence case for psi_fun works!", {
  expect_equal(psi_fun(w_j = c(1, 3, 4),
                       mu = c(4, 13, 45),
                       Sigma =matrix(c(400, 0, 0, 0, 100, 39, 0, 39, 45),
                                     ncol = 3, nrow = 3),
                       params = list(delta_par = c(10, 4, 0.75)),
                       ind = "FALSE"),
               c(1.332620e-12, 2.784265e-13, 4.334743e-12))
})

