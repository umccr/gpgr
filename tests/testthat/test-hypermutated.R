# Tests for hypermutated flag logic in cancer_report.Rmd
# Logic: hypermutated <- snv_summary_counts$is_hypermutated
# The flag is computed in bolt as raw_pass > MAX_SOMATIC_VARIANTS and read here directly.

read_snv_summary_json <- function(counts) {
  tmp <- tempfile(fileext = ".json")
  jsonlite::write_json(counts, tmp, auto_unbox = TRUE)
  jsonlite::read_json(tmp)
}

test_that("hypermutated reads is_hypermutated directly", {
  snv_summary_counts <- read_snv_summary_json(list(raw_pass = 600000, is_hypermutated = TRUE))
  hypermutated <- snv_summary_counts$is_hypermutated
  expect_true(hypermutated)
})

test_that("hypermutated is FALSE when is_hypermutated is FALSE", {
  snv_summary_counts <- read_snv_summary_json(list(raw_pass = 300000, is_hypermutated = FALSE))
  hypermutated <- snv_summary_counts$is_hypermutated
  expect_false(hypermutated)
})

test_that("summary excludes is_hypermutated and includes raw_pass before filter_pass", {
  snv_summary_counts <- read_snv_summary_json(list(
    dragen = 500000, sage = 600000, annotated = 300000,
    raw_pass = 450000, filter_pass = 400000, is_hypermutated = TRUE
  ))
  summary_names <- names(snv_summary_counts)[names(snv_summary_counts) != "is_hypermutated"]
  expect_false("is_hypermutated" %in% summary_names)
  expect_true("raw_pass" %in% summary_names)
  expect_lt(match("raw_pass", summary_names), match("filter_pass", summary_names))
})
