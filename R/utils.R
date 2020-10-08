#' Does the input file look like a VCF?
#'
#' Quickly checks that the input file has a 'vcf' or 'vcf.gz' suffix, and that
#' the column names correspond to typical VCF column headers.
#'
#' @param x Path to file.
#'
#' @return TRUE if the input file is inferred to be a VCF, FALSE otherwise.
#'
#' @examples
#' x <- system.file("extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' (y <- is_vcf(x))
#' tmp_file <- tempfile(pattern = "fakeFile", fileext = "vcf")
#' writeLines(c("col1\tcol2\tcol3", "1\t2\t3"), con = tmp_file)
#' (z <- is_vcf(tmp_file))
#'
#' @testexamples
#' expect_true(y)
#' expect_false(z)
#'
#' @export
is_vcf <- function(x) {
  assertthat::assert_that(length(x) == 1, file.exists(x))
  proper_suffix <- grepl("vcf.gz$", x) | grepl("vcf$", x)
  if (!proper_suffix) {
    message(glue::glue("{x} does not have a vcf or vcf.gz suffix."))
    return(FALSE)
  }

  d <- readr::read_tsv(
    file = x,
    comment = "##",
    col_types = readr::cols(.default = "c"),
    n_max = 1)

  vcf_cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
                "FILTER", "INFO", "FORMAT")
  d_cols <- colnames(d)

  if (!((length(d_cols) >= length(vcf_cols)) &
        (all(d_cols[1:length(vcf_cols)] == vcf_cols)))) {
    message(glue::glue("VCF main column names are incorrect. ",
                       "They are:\n{paste(d_cols, collapse = ', ' )}.\n",
                       "They should include at the beginning:\n",
                       "{paste(vcf_cols, collapse = ', ')}."))
    return(FALSE)
  }
  return(TRUE)
}


#' Does the VCF file contain any variants?
#'
#' Checks if the VCF file contains any variants.
#'
#' @param x Path to VCF file.
#'
#' @return TRUE if there is at least one variant in the VCF, FALSE otherwise.
#'
#' @examples
#' x <- system.file("extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' (y <- vcf_is_empty(x))
#'
#' tmp1 <- tempfile(pattern = "fakeFile", fileext = "vcf")
#' writeLines(c("col1\tcol2\tcol3", "1\t2\t3"), con = tmp1)
#' \dontrun{
#' vcf_is_empty(tmp1)
#' }
#'
#' vcf_cols <- c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
#'               "FILTER", "INFO", "FORMAT")
#' tmp2 <- tempfile(pattern = "fakeFile", fileext = "vcf")
#' writeLines(paste(vcf_cols, collapse = "\t"), con = tmp2)
#' (z <- vcf_is_empty(tmp2))
#'
#' tmp3 <- tempfile(pattern = "fakeFile", fileext = "FOO")
#' writeLines(paste(vcf_cols, collapse = "\t"), con = tmp3)
#' \dontrun{
#' vcf_is_empty(tmp3)
#' }
#'
#' @testexamples
#' expect_false(y)
#' expect_error(vcf_is_empty(tmp1))
#' expect_error(vcf_is_empty(tmp3))
#'
#' @export
vcf_is_empty <- function(x) {
  assertthat::assert_that(is_vcf(x))
  d <- readr::read_tsv(
    file = x,
    comment = "##",
    col_types = readr::cols(.default = "c"),
    n_max = 1)

  if (nrow(d) != 1) {
    message(glue::glue("{x} does not contain any variants."))
    return(TRUE)
  }
  return(FALSE)
}
