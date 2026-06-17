test_that("ncei works", {
  expect_no_error(get_storms("ncei"))
})
test_that("hurdat works", {
  expect_no_error(get_storms("hurdat"))
})
test_that("csv works", {
  expect_no_error(get_storms("E:/OneDrive - USDA/Hurricanes/ibtracs.ALL.list.v04r01.csv"))
})
test_that("nc works", {
  expect_no_error(get_storms("E:/OneDrive - USDA/Hurricanes/IBTrACS.NA.v04r01.nc"))
})

