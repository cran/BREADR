test_that("saveSLICES errors", {
  empty_tibble <- tibble::tribble(
    ~name, ~value
  )
  expect_error(saveSLICES(empty_tibble))
})
