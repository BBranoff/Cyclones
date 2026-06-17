allstorms <- get_storms("E:/OneDrive - USDA/Hurricanes/ibtracs.ALL.list.v04r01.csv")
fewstorms <-  get_storms(ib_filt="last3years")
test_that("build_models works", {
  expect_no_error(build_models(allstorms))
})
test_that("too few build_models works", {
  expect_no_error(build_models(fewstorms))
})
