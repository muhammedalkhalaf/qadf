test_that("qadf works with random walk", {
  set.seed(123)
  y <- cumsum(rnorm(200))
  
  result <- qadf(y)
  
  expect_s3_class(result, "qadf")
  expect_true("results" %in% names(result))
  expect_true("qks" %in% names(result))
  expect_true(is.numeric(result$qks))
  expect_true(result$qks > 0)
})

test_that("qadf works with stationary series", {
  set.seed(456)
  y <- arima.sim(list(ar = 0.5), n = 200)
  
  result <- qadf(y)
  
  expect_s3_class(result, "qadf")
  expect_true(nrow(result$results) == 9)  # Default 9 quantiles
})

test_that("qadf respects model argument", {
  set.seed(789)
  y <- cumsum(rnorm(200))
  
  result_c <- qadf(y, model = "c")
  result_ct <- qadf(y, model = "ct")
  result_nc <- qadf(y, model = "nc")
  
  expect_equal(result_c$model, "c")
  expect_equal(result_ct$model, "ct")
  expect_equal(result_nc$model, "nc")
})

test_that("qadf respects tau argument", {
  set.seed(101)
  y <- cumsum(rnorm(200))
  
  result <- qadf(y, tau = c(0.25, 0.5, 0.75))
  
  expect_equal(nrow(result$results), 3)
  expect_equal(result$results$tau, c(0.25, 0.5, 0.75))
})

test_that("qadf critical values are correct shape", {
  cv <- qadf:::qadf_critical_values(0.5, "c")
  
  expect_length(cv, 3)
  expect_true(all(names(cv) == c("1%", "5%", "10%")))
  expect_true(cv["1%"] < cv["5%"])
  expect_true(cv["5%"] < cv["10%"])
})

test_that("qadf bandwidth is positive", {
  h <- qadf:::qadf_bandwidth(0.5, 200, hs = TRUE)
  
  expect_true(h > 0)
  expect_true(h < 1)
})

test_that("print and summary methods work", {
  set.seed(202)
  y <- cumsum(rnorm(100))
  result <- qadf(y)
  
  expect_output(print(result), "Quantile Augmented Dickey-Fuller")
  expect_output(summary(result), "Full Results")
})

test_that("qadf handles short series with error", {
  y <- rnorm(10)
  
  expect_error(qadf(y), "too short")
})

test_that("qadf_bootstrap works", {
  skip_on_cran()  # Skip on CRAN due to time
  
  set.seed(303)
  y <- cumsum(rnorm(100))
  result <- qadf(y)
  
  result_boot <- qadf_bootstrap(result, nboot = 50, seed = 123)
  
  expect_true("qks_pvalue" %in% names(result_boot))
  expect_true(result_boot$qks_pvalue >= 0)
  expect_true(result_boot$qks_pvalue <= 1)
})
