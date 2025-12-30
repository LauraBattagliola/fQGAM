test_that("prep_fqgam_data creates correct structure", {
  # Create simple test data
  n <- 20
  S <- 10
  y <- rnorm(n)
  id <- rep(1:4, each = 5)
  time <- rep(1:5, 4)
  Xobs <- matrix(rnorm(n * S), n, S)
  
  result <- prep_fqgam_data(y, id, time, Xobs)
  
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), n)
  expect_true("y" %in% names(result))
  expect_true("id" %in% names(result))
  expect_true("time" %in% names(result))
  expect_true("Xobs" %in% names(result))
  expect_true("sMat" %in% names(result))
  expect_true("tMat" %in% names(result))
  expect_s3_class(result$id, "factor")
})

test_that("prep_fqgam_data validates input dimensions", {
  n <- 20
  S <- 10
  y <- rnorm(n)
  id <- rep(1:4, each = 5)
  time <- rep(1:5, 4)
  Xobs <- matrix(rnorm(n * S), n, S)
  
  # Wrong y length

  expect_error(prep_fqgam_data(rnorm(5), id, time, Xobs))
  
  # Wrong id length
  expect_error(prep_fqgam_data(y, 1:5, time, Xobs))
  
  # Wrong time length
  expect_error(prep_fqgam_data(y, id, 1:5, Xobs))
})

test_that("prep_fqgam_data handles custom s_range", {
  n <- 20
  S <- 10
  y <- rnorm(n)
  id <- rep(1:4, each = 5)
  time <- rep(1:5, 4)
  Xobs <- matrix(rnorm(n * S), n, S)
  
  result <- prep_fqgam_data(y, id, time, Xobs, s_range = c(0, 2))
  
  expect_equal(min(result$sMat), 0)
  expect_equal(max(result$sMat), 2)
})

test_that("boot_ci validates input", {
  # Create mock bootstrap result
  B <- 10
  maxT <- 5
  mock_boot <- list(
    predDiff = matrix(rnorm(B * maxT), B, maxT),
    tau = 0.5
  )
  class(mock_boot) <- c("fqgam_boot", "list")
  
  original_pred <- rnorm(maxT)
  
  # Should work
  result <- boot_ci(mock_boot, original_pred)
  expect_s3_class(result, "fqgam_ci")
  
  # Wrong length should error
  expect_error(boot_ci(mock_boot, rnorm(3)))
})

test_that("boot_ci computes correct CI methods", {
  B <- 100
  maxT <- 5
  
  # Create mock bootstrap result with known properties
  set.seed(123)
  mock_boot <- list(
    predDiff = matrix(rnorm(B * maxT), B, maxT),
    tau = 0.5
  )
  class(mock_boot) <- c("fqgam_boot", "list")
  
  original_pred <- rep(0, maxT)
  
  # Normal method
  ci_normal <- boot_ci(mock_boot, original_pred, method = "normal")
  expect_equal(ci_normal$method, "normal")
  expect_true(all(ci_normal$ci_lower < ci_normal$ci_upper))
  
  # Percentile method
  ci_pct <- boot_ci(mock_boot, original_pred, method = "percentile")
  expect_equal(ci_pct$method, "percentile")
  
  # Basic method
  ci_basic <- boot_ci(mock_boot, original_pred, method = "basic")
  expect_equal(ci_basic$method, "basic")
})
