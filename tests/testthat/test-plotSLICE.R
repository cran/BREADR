test_that("plotSLICE errors", {
  empty_tibble <- tibble::tribble(
    ~name, ~value
  )
  expect_error(plotSLICE(empty_tibble))
  false_example1 <- relatedness_example %>% dplyr::rename(bob = pair)
  expect_error(plotSLICE(false_example1))
  expect_error(plotSLICE(relatedness_example))
  expect_error(plotSLICE(relatedness_example, row = 100))
  expect_error(plotSLICE(relatedness_example, row = 1, labels = "a"))
})
