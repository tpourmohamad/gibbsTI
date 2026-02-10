test_that("gibbsTI errors on bad input", {

  expect_error(gibbsTI("a"))
  expect_error(gibbsTI(c(1,2)))
  expect_error(gibbsTI(c(1,2,3,4,5), content = 1.2))

})
