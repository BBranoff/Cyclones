fewstorms <-  get_storms(ib_filt="last3years")
test_that("get_winds works", {
  expect_no_error(lapply(storms[1:3],get_winds,methods="all"))
})
