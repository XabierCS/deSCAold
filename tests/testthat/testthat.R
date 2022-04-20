test_that('test SlideWindowMean',{
  v <-  c(2, 4, 6, 8, 10)
  b <- c(4, 6, 8, 9, 10)
  expect_equal(SlideWindowMean(v,2),b)
})