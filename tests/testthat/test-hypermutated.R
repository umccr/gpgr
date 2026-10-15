# Tests for hypermutated and filtering-type logic in cancer_report.Rmd
# Logic:
#   hypermutated <- snv_summary_counts$is_hypermutated
#   filtering_type <- snv_summary_counts$filtering_type (defaults to "none")
# Both fields are computed in bolt's annotate step and read here directly.

read_snv_summary_json <- function(counts) {
  tmp <- tempfile(fileext = ".json")
  jsonlite::write_json(counts, tmp, auto_unbox = TRUE)
  jsonlite::read_json(tmp)
}

filtering_type_labels <- c(
  "none" = "None",
  "max_variants_non_pass" = "Only passing DRAGEN variants included",
  "max_variants_gnomad" = "Common population variants excluded",
  "max_variants_non_cancer_genes" = "Only cancer gene variants included",
  "unresolved" = "No viable filtering strategy available"
)

filtering_type_display <- function(counts) {
  filtering_type <- counts$filtering_type
  if (is.null(filtering_type)) filtering_type <- "none"
  unname(filtering_type_labels[filtering_type])
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

test_that("summary excludes is_hypermutated and filtering_type, includes raw_pass before filter_pass", {
  snv_summary_counts <- read_snv_summary_json(list(
    dragen = 500000, sage = 600000, annotated = 300000,
    raw_pass = 450000, filter_pass = 400000,
    is_hypermutated = TRUE, filtering_type = "max_variants_gnomad"
  ))
  summary_names <- names(snv_summary_counts)[!names(snv_summary_counts) %in% c("is_hypermutated", "filtering_type")]
  expect_false("is_hypermutated" %in% summary_names)
  expect_false("filtering_type" %in% summary_names)
  expect_true("raw_pass" %in% summary_names)
  expect_lt(match("raw_pass", summary_names), match("filter_pass", summary_names))
})

test_that("filtering_type maps every machine value to its display label", {
  expect_equal(filtering_type_display(list(filtering_type = "none")), "None")
  expect_equal(filtering_type_display(list(filtering_type = "max_variants_non_pass")), "Only passing DRAGEN variants included")
  expect_equal(filtering_type_display(list(filtering_type = "max_variants_gnomad")), "Common population variants excluded")
  expect_equal(filtering_type_display(list(filtering_type = "max_variants_non_cancer_genes")), "Only cancer gene variants included")
  expect_equal(filtering_type_display(list(filtering_type = "unresolved")), "No viable filtering strategy available")
})

test_that("filtering_type defaults to none when absent", {
  expect_equal(filtering_type_display(list()), "None")
})
