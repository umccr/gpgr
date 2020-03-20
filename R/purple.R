#' Read PURPLE CNV Gene File
#'
#' Reads the `purple.cnv.gene.tsv` file, which summarises copy number
#' alterations of each gene in the HMF panel
#' (see [this table](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator#gene-copy-number-file)).
#'
#' @param x Path to `purple.cnv.gene.tsv` file.
#' @param v PURPLE version (default: 2.39). Used to determine the column names.
#'
#' @return The input file in a tibble.
#'
#'
#' @examples
#' x <- system.file("extdata/purple/v2.39/purple.cnv.gene.tsv", package = "gpgr")
#' p <- read_purple_cnv_gene(x, "2.39")
#' p
#'
#'
#' @testexamples
#' expect_equal(colnames(p)[ncol(p)], "minMinorAllelePloidy")
#' expect_error(read_purple_cnv_gene(x, "2.38"))
#'
#' @export
read_purple_cnv_gene <- function(x, v = "2.39") {

  # Use utils::compareVersion(a, b) which gives -1/0/1, if a < b, a == b, a > b
  config <- function(v) {
    current <- list(
      nm = c("chromosome", "start", "end", "gene", "minCopyNumber", "maxCopyNumber",
             "unused", "somaticRegions", "germlineHomDeletionRegions", "germlineHetToHomDeletionRegions",
             "transcriptId", "transcriptVersion", "chromosomeBand", "minRegions",
             "minRegionStart", "minRegionEnd", "minRegionStartSupport", "minRegionEndSupport",
             "minRegionMethod", "minMinorAllelePloidy"),
      type = "ciicddcdddcccdiicccd")
    l <- list(
      "2.39" = current,
      "2.45" = list(nm = c(current$nm, "foo", "bar"), type = paste0(current$type, "cc"))
    )
    if (utils::compareVersion("2.45", v) <= 0) { # 2.45, 2.46...
      return(l[["2.45"]])
    } else if (utils::compareVersion("2.39", v) <= 0) { # 2.39, 2.40... 2.44
      return(l[["2.39"]])
    } else {
      stop(glue::glue("You've specified PURPLE version {v}, which is older than 2.39.",
                      "Please update to newer PURPLE version."))
    }
  }
  conf <- config(v)
  d <- readr::read_tsv(x, col_types = conf$type)
  assertthat::assert_that(all(colnames(d) == conf$nm))
  d
}
