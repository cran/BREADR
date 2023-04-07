test_that("test_degree errors", {
  empty_tibble <- tibble::tribble(
    ~name, ~value
  )
  false_example <- relatedness_example %>% dplyr::rename(bob = pair)
  expect_error(test_degree(empty_tibble))
  expect_error(test_degree(false_example))
  expect_error(test_degree(relatedness_example, row = 100))
  expect_error(test_degree(relatedness_example, row = 1, degree = "bob"))
  expect_equal(test_degree(relatedness_example, row = 1, degree = 1), 2.667e-5,
               tolerance = 0.0001)
})
