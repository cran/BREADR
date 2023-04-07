test_that("split_line works", {
  expect_error(split_line("1\t1"))
  expect_equal(as.character(split_line("rs3094315\t1\t0.0\t752566\tG\tA")[1,1]), "rs3094315")
  expect_equal(as.numeric(split_line("rs3094315\t1\t0.0\t752566\tG\tA")[1,2]), 1)
})
