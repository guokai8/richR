test_that("reports prefer native Quartz even when Cairo is advertised", {
  skip_if_not(capabilities("png"), "the png device is unavailable")
  native_type <- if (capabilities("aqua")) "quartz" else {
    if (capabilities("cairo")) "cairo" else {
      if (.Platform$OS.type == "windows") "windows" else "Xlib"
    }
  }
  png <- grDevices::png
  testthat::local_mocked_bindings(
    capabilities = function(what) stats::setNames(rep(TRUE, length(what)), what),
    .package = "base")
  testthat::local_mocked_bindings(
    png = function(..., type) {
      expect_identical(type, "quartz")
      png(..., type = native_type)
    }, .package = "grDevices")
  before <- grDevices::dev.cur()
  drawn <- richR:::.report_draw(.rp_fq()$q1, FALSE, 30)
  expect_identical(grDevices::dev.cur(), before)
  bytes <- jsonlite::base64_dec(drawn$png)
  expect_gt(length(bytes), 100L)
  expect_identical(bytes[1:8], as.raw(c(137, 80, 78, 71, 13, 10, 26, 10)))
})
