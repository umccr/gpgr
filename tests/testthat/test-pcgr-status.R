# Tests for pcgr_format_categories() (sash #52).
# Previously this logic lived inline in cancer_report.Rmd and was
# reimplemented here just to get coverage. It's now extracted into
# R/umccrise.R and exported by gpgr, so these tests exercise the real
# function directly (see also the roxytest-generated
# test-roxytest-testexamples-umccrise.R, which covers its @examples).

test_that("pcgr_format_categories parses a single category", {
  out <- pcgr_format_categories("N|intronic|difficult:80321")
  expect_equal(out, "- non-coding / intronic / difficult region: 80,321 variants")
})

test_that("pcgr_format_categories parses multiple comma-separated categories", {
  out <- pcgr_format_categories("N|intronic|difficult:80321,2|impacts_other|giab_conf:5000")
  expect_equal(
    out,
    paste(
      "- non-coding / intronic / difficult region: 80,321 variants",
      "- tier 2 / other VEP consequence / GIAB confident region: 5,000 variants",
      sep = "\n"
    )
  )
})

test_that("pcgr_format_categories joins multi-element vector input with commas, not fused together", {
  # bolt may emit this as a character vector rather than one pre-joined string
  out <- pcgr_format_categories(c("N|intronic|difficult:80321", "2|impacts_other|giab_conf:5000"))
  expect_equal(
    out,
    paste(
      "- non-coding / intronic / difficult region: 80,321 variants",
      "- tier 2 / other VEP consequence / GIAB confident region: 5,000 variants",
      sep = "\n"
    )
  )
})

test_that("pcgr_format_categories returns empty string for empty input", {
  expect_equal(pcgr_format_categories(""), "")
  expect_equal(pcgr_format_categories(character(0)), "")
})

test_that("pcgr_format_categories skips a malformed entry instead of erroring", {
  # missing the region:count segment on the second entry
  out <- pcgr_format_categories("N|intronic|difficult:80321,malformed_entry")
  expect_equal(out, "- non-coding / intronic / difficult region: 80,321 variants")
})

test_that("pcgr_format_categories skips an entry with a non-numeric count instead of erroring", {
  out <- pcgr_format_categories("N|intronic|difficult:not_a_number")
  expect_equal(out, "")
})

test_that("pcgr_format_categories falls back to the raw code for unknown labels", {
  out <- pcgr_format_categories("9|unknown_impact|unknown_region:100")
  expect_equal(out, "- 9 / unknown_impact / unknown_region: 100 variants")
})
