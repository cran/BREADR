test_that("plotLOAF errors", {
  empty_tibble <- tibble::tribble(
    ~name, ~value
  )
  expect_error(plotLOAF(empty_tibble))
  false_example1 <- relatedness_example %>% dplyr::rename(bob = pair)
  expect_error(plotLOAF(false_example1))
  expect_error(plotLOAF(kin_example, nsnps_cutoff = -1))
})
