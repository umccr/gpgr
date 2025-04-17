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

  readr::local_edition(1)
  d <- readr::read_tsv(
    file = x,
    comment = "##",
    col_types = readr::cols(.default = "c"),
    n_max = 1
  )

  # COMMIT NOTE: FORMAT FIELD IS NOT MANDITORY

  vcf_cols <- c(
    "#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
    "FILTER", "INFO"
  )
  d_cols <- colnames(d)

  if (!((length(d_cols) >= length(vcf_cols)) &
    (all(d_cols[1:length(vcf_cols)] == vcf_cols)))) {
    message(glue::glue(
      "VCF main column names are incorrect. ",
      "They are:\n{paste(d_cols, collapse = ', ' )}.\n",
      "They should include at the beginning:\n",
      "{paste(vcf_cols, collapse = ', ')}."
    ))
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
#' vcf_cols <- c(
#'   "#CHROM", "POS", "ID", "REF", "ALT", "QUAL",
#'   "FILTER", "INFO", "FORMAT"
#' )
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
#' @testexamples
#' expect_false(a)
#' expect_true(b)
#' expect_true(c)
#'
#' @export
tsv_is_empty <- function(x, comment = "##", col_types = readr::cols(.default = "c"), n_max = 1, ...) {
  readr::local_edition(1)
  d <- readr::read_tsv(
    file = x,
    comment = comment,
    col_types = col_types,
    n_max = n_max,
    ...
  )

  if (nrow(d) != n_max) {
    message(glue::glue("{x} does not contain any non-header rows."))
    return(TRUE)
  }
  return(FALSE)
}

#' Write a data frame to a tab delimited gzipped file
#'
#' Writes a data frame to a tab delimited gzipped file.
#'
#' @param x A data frame or tibble to write to disk.
#' @param file File or connection to write to (should end in '.gz').
#' @param ... Additional arguments passed to [readr::write_tsv()].
#'
#' @return Returns the input `x` invisibly.
#'
#' @export
write_tsvgz <- function(x, file, ...) {
  assertthat::assert_that(endsWith(file, ".gz"), inherits(x, "data.frame"))
  mkdir(dirname(file))
  readr::write_tsv(x = x, file = file, ...)
}

#' Write gzipped JSON
#'
#' Serializes an object to JSON and writes to a gzipped file.
#'
#' @param x An object to be serialized to JSON.
#' @param path File on disk (should end in '.gz').
#' @param ... Additional arguments passed to [jsonlite::write_json()]
#'
#' @export
write_jsongz <- function(x, path, ...) {
  assertthat::assert_that(endsWith(path, ".gz"))
  mkdir(dirname(path))
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
  write_jsongz(x, file.path(maindir, "json", json_path))
}

#' Does R Package Exist
#'
#' Checks if the specified R package exists on the local system.
#'
#' @param p The R package to check for.
#' @return TRUE if package exists, FALSE otherwise.
#'
#' @export
pkg_exists <- function(p) {
  assertthat::assert_that(is.character(p))
  nzchar(system.file(package = p))
}

#' Print current timestamp for logging
#'
#' @return Current timestamp as character.
#' @export
date_log <- function() {
  as.character(paste0("[", as.POSIXct(Sys.time()), "]"))
}

#' Create directory
#'
#' @param d Directory to create.
#'
#' @export
mkdir <- function(d) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
}

# copy recursively
cpdir <- function(from, to) {
  mkdir(to)
  file.copy(from = from, to = to, recursive = TRUE)
}

# https://stackoverflow.com/a/61647053/2169986
mixedrank <- function(x) {
  if (length(x) == 0) {
    return(x)
  }
  order(gtools::mixedorder(x))
}

is_url <- function(x) {
  grepl("(http|https)://[a-zA-Z0-9./?=_%:-]*", x)
}

bcftools_installed <- function() {
  system("bcftools -v", ignore.stdout = TRUE) == 0
}

#' Parse VCF with bcftools
#'
#' Parse VCF with bcftools.
#' Uses bcftools under the hood to do the heavy lifting with field splitting,
#' then converts the parsed character vector to a tibble.
#'
#' For VCFs with 0 variants, returns a tibble with 0 rows and proper number of
#' columns.
#'
#' @param vcf VCF with one or more samples.
#' @param only_pass Keep PASS variants only (def: TRUE).
#' @param alias Substitute sample names with S1/S2/... alias (def: TRUE).
#'
#' @return A tibble with all the main, FORMAT, and INFO fields detected in
#' the VCF header as columns.
#' @export
bcftools_parse_vcf <- function(vcf, only_pass = TRUE, alias = TRUE) {
  assertthat::assert_that(is.logical(only_pass), length(only_pass) == 1)
  assertthat::assert_that(
    bcftools_installed(),
    msg = "bcftools needs to be on the PATH."
  )

  if (is_url(vcf)) {
    vcf <- glue::glue("'{vcf}'")
  }

  # 1) grab the header
  cmd_header <- glue::glue("bcftools view -h {vcf}")
  h <- system(cmd_header, intern = TRUE)

  # 2) main fields + sample aliases
  main <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER")
  get_samples <- function() {
    hdr_last <- strsplit(h[length(h)], "\t")[[1]]
    # drop CHROM–FILTER + INFO + FORMAT
    samples <- hdr_last[-seq_len(length(main) + 2)]
    aliases <- if (alias) paste0("S", seq_along(samples)) else samples
    list(samples = samples, n = length(samples), aliases = aliases)
  }
  samp <- get_samples()

  # 3) function to split INFO or FORMAT lines
  split_hdr <- function(pat) {
    h[grepl(pat, h, fixed = TRUE)] |>
      sub(paste0("^", pat), "", _)          |>  # strip the prefix
      tibble::as_tibble_col(column_name = "x") |>
      tidyr::separate_wider_delim(
        "x",
        delim = ",",
        names = c("ID", "Number", "Type", "Description"),
        too_many = "merge"
      ) |>
      dplyr::mutate(
        ID          = sub("^ID=",          "", .data$ID),
        Number      = sub("^Number=",      "", .data$Number),
        Type        = sub("^Type=",        "", .data$Type),
        Description = sub('^Description="(.*)">?$', "\\1", .data$Description)
      )
  }

  fmt <- split_hdr("##FORMAT=<") |> dplyr::pull(.data$ID)
  info <- split_hdr("##INFO=<")   |> dplyr::pull(.data$ID)

  # 4) build the bcftools query string
  main_cols  <- paste0("%", main, collapse = "\\t")
  info_cols  <- if (length(info) == 0) "" else paste0("%INFO/", info, collapse = "\\t")
  fmt_cols   <- paste0("[\\t", paste0("%", fmt, collapse = "\\t"), "]\\n")
  q          <- paste0(main_cols, "\\t", info_cols, fmt_cols)
  include_pass <- if (only_pass) "-i 'FILTER=\"PASS\" || FILTER=\".\"'" else ""
  cmd_body     <- glue::glue("bcftools query -f \"{q}\" {vcf} {include_pass}")

  # 5) build the vector of desired names
  desired_names <- c(
    main,
    if (length(info) > 0) paste0("INFO_", info),
    paste0(rep(samp$aliases, each = length(fmt)), "_", fmt)
  )

  # 6) preview to get actual column count
  suppressWarnings({
    preview <- data.table::fread(
      cmd        = cmd_body,
      sep        = "\t",
      na.strings = ".",
      nrows      = 1,
      header     = FALSE
    )
  })
  actual_ncol <- ncol(preview)

  # handle empty VCF
  if (nrow(preview) == 0) {
    return(empty_tbl(cnames = desired_names))
  }

  # ensure we have at least actual_ncol names
  if (length(desired_names) < actual_ncol) {
    stop(
      sprintf(
        "Expected at least %d column names, but only have %d. Check your header fields.",
        actual_ncol, length(desired_names)
      )
    )
  }

  # 7) read full table and assign exactly the right number of names
  df <- data.table::fread(
    cmd        = cmd_body,
    sep        = "\t",
    header     = FALSE,
    na.strings = ".",
    data.table = FALSE
  )
  colnames(df) <- desired_names[seq_len(actual_ncol)]
  df |> tibble::as_tibble()
}
