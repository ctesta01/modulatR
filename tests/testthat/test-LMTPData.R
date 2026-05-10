library(testthat)

# tests for the LMTPData class

# test that the constructor works

testthat::test_that("LMTPData constructor works", {


  expect_no_error({
    n0 <- 10

    ds <- LMTPData$new(
      data = data.frame(
        L = rnorm(n = n0),
        A = rnorm(n = n0),
        Y = rnorm(n = n0)),
      L_cols = 'L',
      A_cols = 'A',
      Y_col = 'Y'
    )
  })

  expect_no_error({

    n0 <- 10

    ds <- LMTPData$new(
      data = tibble::tibble(
        W0 = rnorm(n0),
        W1 = rnorm(n0),
        Lt1_1 = rnorm(n0),
        Lt1_2 = rpois(n0, 5),
        At1 = rbinom(n0, prob = plogis(Lt1_1 + Lt1_2 + W0 + W1), size = 1),
        Lt2_1 = rnorm(n0, mean = Lt1_1),
        Lt2_2 = rpois(n0, lambda = Lt1_2 + 1),
        At2 = rbinom(n0, prob = plogis(Lt2_1 + Lt2_2 + At1), size = 1),
        Y = W0 + W1 + Lt1_1 + Lt1_2 + Lt2_1 + Lt2_2 + rnorm(n0)
      ),
      W_cols = c('W0', 'W1'),
      L_cols = list(t1 = c('Lt1_1', 'Lt1_2'), t2 = c('Lt2_1', 'Lt2_2')),
      A_cols = c('At1', 'At2'),
      Y_col = 'Y')
  })



})
