#' Read VIRUSBreakend Summary File
#'
#' Reads the `virusbreakend.vcf.summary.tsv` file.
#'
#' @param x Path to `virusbreakend.vcf.summary.tsv` file.
#'
#' @return List with two elements:
#' * `tab`: Tibble containing data.
#' * `descr`: Description of tibble columns.
#'
#' @examples
#' x <- system.file("extdata/virusbreakend/virusbreakend.vcf.summary.tsv", package = "gpgr")
#' (vb <- virusbreakend_summary_read(x))
#' @testexamples
#' expect_equal(colnames(vb$tab)[ncol(vb$tab)], "QC")
#'
#' @export
virusbreakend_summary_read <- function(x) {
  nm <- c(
    "taxid_genus" = "c",
    "name_genus" = "c",
    "reads_genus_tree" = "i",
    "taxid_species" = "c",
    "name_species" = "c",
    "reads_species_tree" = "i",
    "taxid_assigned" = "c",
    "name_assigned" = "c",
    "reads_assigned_tree" = "i",
    "reads_assigned_direct" = "i",
    "reference" = "c",
    "reference_taxid" = "c",
    "reference_kmer_count" = "i",
    "alternate_kmer_count" = "i",
    "rname" = "c",
    "startpos" = "i",
    "endpos" = "i",
    "numreads" = "i",
    "covbases" = "d",
    "coverage" = "d",
    "meandepth" = "d",
    "meanbaseq" = "d",
    "meanmapq" = "d",
    "integrations" = "i",
    "QCStatus" = "c"
  )
  ctypes <- paste(nm, collapse = "")
  virusbreakend_summary <- readr::read_tsv(x, col_types = ctypes)

  if (nrow(virusbreakend_summary) > 0) {
    assertthat::assert_that(ncol(virusbreakend_summary) == length(nm))
    assertthat::assert_that(all(colnames(virusbreakend_summary) == names(nm)))

    virusbreakend_summary <- virusbreakend_summary |>
      dplyr::select(
        Virus = "name_assigned",
        Length = "endpos",
        Reads = "numreads",
        Coverage = "coverage",
        `Mean depth` = "meandepth",
        Intergrations = "integrations",
        QC = "QCStatus",
      )
  }

  descr <- dplyr::tribble(
    ~Column,
    ~Description,
    "Virus",
    "Assigned NCBI taxonomy name of viral reference",
    "Length",
    "Length of viral contig",
    "Reads",
    "Number of reads mapped to adjusted viral reference",
    "Coverage",
    "Percentage of viral positions with at least one read mapped",
    "Mean depth",
    "Mean alignment depth",
    "Integrations",
    "Number of detected integration breakpoints",
    "QC",
    "QC status of viral intergrations",
  )

  list(
    tab = virusbreakend_summary,
    descr = descr
  )
}

#' Read VIRUSBreakend VCF File
#'
#' Reads the `virusbreakend.vcf` file and selects data to present.
#'
#' @param x Path to `virusbreakend.vcf` file.
#'
#' @return List with two elements:
#' * `tab`: Tibble containing selected data.
#' * `descr`: Description of tibble columns.
#'
#' @examples
#' x <- system.file("extdata/virusbreakend/virusbreakend.vcf", package = "gpgr")
#' (vb <- virusbreakend_vcf_read(x))
#' @testexamples
#' expect_equal(colnames(vb$tab)[ncol(vb$tab)], "QC")
#'
#' @export
virusbreakend_vcf_read <- function(x) {
  d <- bedr::read.vcf(x, split.info = TRUE, verbose = FALSE)

  if (nrow(d$vcf) > 0) {
    virusbreakend_integrations <- tibble::as_tibble(d$vcf) |>
      dplyr::select(
        Contig = "CHROM",
        Position = "POS",
        "Fragment support" = "BVF",
        "Fragment support (unmapped)" = "BUM",
        "Softclip read support" = "BSC",
        Reference = "REF",
        Alt = "ALT",
        `Breakend ID` = "ID",
        `Mate ID` = "MATEID",
        QC = "FILTER",
      )
  } else {
    virusbreakend_integrations <- tibble::tibble()
  }

  descr <- dplyr::tribble(
    ~Column,
    ~Description,
    "Contig",
    "Name of contig",
    "Position",
    "Position of breakend in contig",
    "Breakend ID",
    "ID of integration breakend",
    "Mate ID",
    "ID of integration breakend mate",
    "Reference",
    "Reference allele",
    "Alt",
    "Alternative allele",
    "QC",
    "VCF filter values",
    "Fragment support",
    "Total number of fragments supporting breakend",
    "Fragment support (unmapped)",
    "Number of fragments supporting breakend that have one read unmapped",
    "Softclip read support",
    "Number of softclipped reads supporting breakend"
  )

  list(
    tab = virusbreakend_integrations,
    description = descr
  )
}
