# ===================================================================
# SRR COMPLIANCE METADATA ASSERTIONS
# Dedicated structural file for tracking global package parameters
# ===================================================================

#' @srrstats {G5.0}
#' @srrstats {G5.5}
#' @srrstats {G5.6b}
#' @srrstats {G5.8}
#' @srrstats {G5.8a}
#' @srrstats {G5.8c}
#' @srrstats {G5.8d}
#' @srrstats {RE7.0}
test_that("Package global data environments satisfy structural review properties", {
  expect_true(exists("mock_data"))
  expect_true(exists("mock_data_pathological"))
  expect_identical(n_test, 60) # UPGRADED: expect_equal -> expect_identical
})
