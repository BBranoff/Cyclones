fewstorms <-  get_storms(ib_filt="last3years")
test_that("make_extents works", {
  expect_no_error(lapply(storms[1:3],make_extents))
})
