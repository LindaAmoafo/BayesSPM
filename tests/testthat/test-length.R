#context("checking length")

test_that("checking length of group with Y", {
  Y = mvtnorm::rmvnorm(40,mean=rep(0,50))
  grp=rep(c(1,2,3,4),each=10)
  time=1:50
  out <- SPMBayes(Y,group=grp, time, threshold=0.95, Q.threshold=0.05, Hypothesis="alt")
  expect_equal(ncol(out$SUMMARY), 18)

})

## add a test for ppq in (0,1)
