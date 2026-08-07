# Tests for hrd_results_tabs() (R/umccrise.R), covering both the DRAGEN-present
# and DRAGEN-absent (dragen_res = NULL, e.g. OA-only mode) code paths.
# Column counts below are chosen to satisfy the function's row-padding
# contract (chord/hrdetect tibbles must end up with equal row counts for
# dplyr::bind_cols() to succeed) — not a claim about the real HRDetect shape.

mock_chord_res <- function() {
  list(
    prediction = tibble::tibble(
      sample = "S1", p_hrd = 0.8, p_BRCA1 = 0.1, p_BRCA2 = 0.7,
      hr_status = "HR_deficient", hrd_type = "BRCA2_type",
      remarks_hr_status = "", remarks_hrd_type = ""
    )
  )
}

mock_hrdetect_res <- function() {
  tibble::tibble(
    sample = "S1", del.mh.prop = 0.1, SNV3 = 0.2, SV3 = 0.1, SV5 = 0.05,
    hrdetect_prob = 0.9, Probability = 0.9, intercept = 0.1, dummy = 0.1
  )
}

mock_dragen_res <- function() {
  res <- gpgr::dragen_hrd(NULL)
  res$HRD <- "0.5"
  res$LOH <- "10"
  res$TAI <- "5"
  res$LST <- "3"
  res
}

test_that("hrd_results_tabs omits the DRAGEN column when dragen_res is NULL", {
  res <- hrd_results_tabs(mock_hrdetect_res(), mock_chord_res(), dragen_res = NULL)
  expect_false("DRAGEN" %in% colnames(res$hrd_results_tab))
  expect_false("results_dragen" %in% colnames(res$hrd_results_tab))
  expect_true(all(c("CHORD", "HRDetect") %in% colnames(res$hrd_results_tab)))
  expect_s3_class(res$hrd_results_gt, "gt_tbl")
})

test_that("hrd_results_tabs includes the DRAGEN column when dragen_res is supplied", {
  res <- hrd_results_tabs(mock_hrdetect_res(), mock_chord_res(), dragen_res = mock_dragen_res())
  expect_true(all(c("DRAGEN", "results_dragen", "CHORD", "HRDetect") %in% colnames(res$hrd_results_tab)))
  expect_s3_class(res$hrd_results_gt, "gt_tbl")
})

test_that("hrd_results_tabs renders to HTML without error in both branches", {
  no_dragen <- hrd_results_tabs(mock_hrdetect_res(), mock_chord_res(), dragen_res = NULL)
  with_dragen <- hrd_results_tabs(mock_hrdetect_res(), mock_chord_res(), dragen_res = mock_dragen_res())
  expect_no_error(gt::as_raw_html(no_dragen$hrd_results_gt))
  expect_no_error(gt::as_raw_html(with_dragen$hrd_results_gt))
})
