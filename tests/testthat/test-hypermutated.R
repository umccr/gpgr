# Tests for hypermutated flag logic in cancer_report.Rmd
# Logic: hypermutated <- snv_summary_counts$sage != snv_summary_counts$annotated
# Hypermutated is inferred indirectly: if sage != annotated, bolt's hypermutated
# logic ran and applied variant selection.

read_snv_summary_json <- function(counts) {
  tmp <- tempfile(fileext = ".json")
  jsonlite::write_json(counts, tmp, auto_unbox = TRUE)
  jsonlite::read_json(tmp)
}

test_that("hypermutated is TRUE when sage != annotated", {
  snv_summary_counts <- read_snv_summary_json(list(sage = 600000, annotated = 300000))
  hypermutated <- snv_summary_counts$sage != snv_summary_counts$annotated
  expect_true(hypermutated)
})

test_that("hypermutated is FALSE when sage == annotated", {
  snv_summary_counts <- read_snv_summary_json(list(sage = 300000, annotated = 300000))
  hypermutated <- snv_summary_counts$sage != snv_summary_counts$annotated
  expect_false(hypermutated)
})

test_that("hypermutated is TRUE when annotated < sage (bolt filtered variants)", {
  snv_summary_counts <- read_snv_summary_json(list(sage = 800000, annotated = 50000))
  hypermutated <- snv_summary_counts$sage != snv_summary_counts$annotated
  expect_true(hypermutated)
})
