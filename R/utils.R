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
  tsv_is_empty(x)
}

#' Does the TSV file contain any non-header rows?
#'
#' Checks if the TSV file contains any non-header rows.
#'
#' @param x Path to TSV file.
#' @param comment Taken from [readr::read_tsv()].
#' @param col_types Taken from [readr::read_tsv()].
#' @param n_max Taken from [readr::read_tsv()].
#' @param ... Other arguments to be passed to [readr::read_tsv()].
#'
#' @return TRUE if there is at least one non-header row in the TSV, FALSE otherwise.
#'
#' @examples
#'
#' tmp1 <- tempfile()
#' writeLines(c("col1\tcol2\tcol3", "1\t2\t3"), con = tmp1)
#' (a <- tsv_is_empty(tmp1))
#'
#' tmp2 <- tempfile()
#' writeLines(c("col1\tcol2\tcol3"), con = tmp2)
#' (b <- tsv_is_empty(tmp2))
#'
#' tmp3 <- tempfile()
#' writeLines(c("##meta1", "##meta2", "col1\tcol2\tcol3"), con = tmp3)
#' (c <- tsv_is_empty(tmp3))
#'
#' @testexamples
#' expect_false(a)
#' expect_true(b)
#' expect_true(c)
#'
#' @export
tsv_is_empty <- function(x, comment = "##", col_types = readr::cols(.default = "c"), n_max = 1, ...) {
  d <- readr::read_tsv(
    file = x,
    comment = comment,
    col_types = col_types,
    n_max = n_max,
    ...)

  if (nrow(d) != n_max) {
    message(glue::glue("{x} does not contain any non-header rows."))
    return(TRUE)
  }
  return(FALSE)
}

#' Make Directory
#'
#' Creates a directory.
#'
#' @param d Directory to create.
#'
#' @return If directory exists, do nothing. If it doesn't, create it and return
#' invisibly a logical indicating if the operation succeeded.
#'
#' @export
mkdir <- function(d) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

#' Write a data frame to a tab delimited gzipped file
#'
#' Writes a data frame to a tab delimited gzipped file.
#'
#' @param x A data frame or tibble to write to disk.
#' @param file File or connection to write to.
#' @param ... Additional arguments passed to [readr::write_tsv()].
#'
#' @return Returns the input `x` invisibly.
#'
#' @export
write_tsvgz <- function(x, file, ...) {
  assertthat::assert_that(endsWith(file, ".gz"), inherits(x, "data.frame"))
  mkdir(dirname(path))
  readr::write_tsv(x = x, file = file, ...)
}

#' Write gzipped JSON
#'
#' Serializes an object to JSON and writes to a gzipped file.
#'
#' @param x An object to be serialized to JSON.
#' @param path File on disk.
#' @param ... Additional arguments passed to [jsonlite::write_json()]
#'
#' @export
write_jsongz <- function(x, path, ...) {
  assertthat::assert_that(endsWith(path, ".gz"))
  mkdir(dirname(json_path))
  gz <- gzfile(path, open = "w")
  jsonlite::write_json(x = x, path = gz, ...)
  close(gz)
}

#' Write a data frame to gzipped TSV and JSON files
#'
#' Writes a data frame to gzipped TSV and JSON files. Files
#' will be written to `maindir/<path>.tsv.gz` and `maindir/json/<path>.json.gz`.
#'
#' @param x The data frame to write.
#' @param path Relative path to write the files to, sans the file extensions -
#' these will be appended appropriately.
#' @param maindir Main directory to write the files to.
#'
#' @export
write_tsvjsongz <- function(x, path, maindir) {
  json_path <- paste0(path, ".json.gz")
  tsv_path <- paste0(path, ".tsv.gz")
  write_tsvgz(x, file.path(maindir, tsv_path))
  write_jsongz(x, file.path(maindir, json_path))
}
