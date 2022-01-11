#' Generate LINX R Markdown HTML Report
#'
#' Generates a LINX R Markdown HTML report. It does so by:
#' 1. recursively copy the table and plot directories into tmp/plot and tmp/table
#' 2. copy the rmd into tmp/linx.Rmd
#' 3. render the rmd inside tmp/
#' 4. return the path to the output HTML
#'
#' @param sample Name of sample.
#' @param table_dir Path to LINX table directory.
#' @param plot_dir Path to LINX plot directory.
#' @param out_file Path to output HTML file (needs '.html' suffix).
#' @param quiet Suppress printing during rendering.
#'
#' @return Path to rendered HTML report.
#' @export
linx_rmd <- function(sample, table_dir, plot_dir, out_file = NULL, quiet = FALSE) {
  assertthat::assert_that(
    is.character(sample), nchar(sample) > 0,
    dir.exists(table_dir), dir.exists(plot_dir),
    quiet %in% c(FALSE, TRUE)
  )
  if (!is.null(out_file)) {
    assertthat::assert_that(
      is.character(out_file),
      tools::file_ext(out_file) == "html"
    )
  } else {
    out_file <- glue::glue("linx_{sample}.html")
  }
  tmp_dir <- tempdir()
  rmd_dir <- system.file("rmd/linx", package = "gpgr")
  cpdir(table_dir, tmp_dir)
  cpdir(plot_dir, tmp_dir)
  cpdir(rmd_dir, tmp_dir)
  rmd_file <- file.path(tmp_dir, "linx", "linx.Rmd")
  out_dir <- dirname(out_file)
  mkdir(out_dir)
  pars <- list(
    table_dir = table_dir,
    plot_dir = plot_dir,
    sample = sample
  )
  # suppress DT large size warning
  options(DT.warn.size = FALSE)
  rmarkdown::render(
    input = rmd_file,
    params = pars,
    output_dir = out_dir,
    output_file = I(out_file),
    run_pandoc = TRUE,
    quiet = quiet
  )
}
