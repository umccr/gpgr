#' Read PURPLE CNV Gene File
#'
#' Reads the `purple.cnv.gene.tsv` file, which summarises copy number
#' alterations of each gene in the HMF panel
#' (see [this table](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator#gene-copy-number-file)).
#'
#' @param x Path to `purple.cnv.gene.tsv` file.
#' @param v PURPLE version (default: 2.39). Used to determine the column names.
#'
#' @return The input file as a tibble.
#'
#' @examples
#' x <- system.file("extdata/purple/v2.39/purple.cnv.gene.tsv", package = "gpgr")
#' p <- read_purple_cnv_gene(x, "2.39")
#' p
#'
#' @testexamples
#' expect_equal(colnames(p)[ncol(p)], "minMinorAllelePloidy")
#' expect_error(read_purple_cnv_gene(x, "2.38"))
#'
#' @export
read_purple_cnv_gene <- function(x, v = "2.39") {

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
    if (numeric_version(v) >= numeric_version("2.45")) {
      return(l[["2.45"]])
    } else if (numeric_version(v) >= numeric_version("2.39")) {
      return(l[["2.39"]])
    } else {
      stop(glue::glue("You've specified PURPLE version {v}, which is older than 2.39.",
                      "Please update to newer PURPLE version."))
    }
  }
  conf <- config(v)
  purple_cnv_gene <- readr::read_tsv(x, col_types = conf$type)
  assertthat::assert_that(ncol(purple_cnv_gene) == length(conf$nm))
  assertthat::assert_that(all(colnames(purple_cnv_gene) == conf$nm))
  purple_cnv_gene
}


#' Process PURPLE CNV Gene File for UMCCRISE
#'
#' Processes the `purple.cnv.gene.tsv` file. Keeps genes that are in the
#' [UMCCR cancer gene list](https://github.com/umccr/genes/blob/893a655801ce92715f05517b5052e4e81904e870/panels/umccr_2019-03-20.tsv)
#' and selects columns of interest.
#'
#' @param x Path to `purple.cnv.gene.tsv` file.
#' @param g Path to gene file containing at least three columns:
#' * `symbol`: gene name (character).
#' * `tumorsuppressor`: is this gene a tumor suppressor (TRUE/FALSE).
#' * `oncogene`: is this gene an oncogene (TRUE/FALSE).
#' @param v PURPLE version (default: 2.39). Used to determine the column names.
#'
#' @return Tibble filtered to genes found in  `g`.
#'
#' @examples
#' x <- system.file("extdata/purple/v2.39/purple.cnv.gene.tsv", package = "gpgr")
#' g <- system.file("extdata/ref/umccr_cancer_genes_2019-03-20.tsv", package = "gpgr")
#' pp <- process_purple_cnv_gene(x, g)
#'
#' @testexamples
#' expect_equal(colnames(pp)[ncol(pp)], "minRegSupportStartEndMethod")
#' expect_error(process_purple_cnv_gene(x, "2.38"))
#'
#' @export
process_purple_cnv_gene <- function(x, g = NULL, v = "2.39") {
  purple_cnv_gene <- read_purple_cnv_gene(x, v)
  if (is.null(g)) {
    g <- system.file("extdata/ref/umccr_cancer_genes_2019-03-20.tsv", package = "gpgr")
  }
  genes <- readr::read_tsv(g, col_types = readr::cols(symbol = "c", oncogene = "l", tumorsuppressor = "l")) %>%
    dplyr::select(.data$symbol, .data$oncogene, .data$tumorsuppressor)
  oncogenes <- genes %>% dplyr::filter(.data$oncogene) %>% dplyr::pull(.data$symbol)
  tsgenes <- genes %>% dplyr::filter(.data$tumorsuppressor) %>% dplyr::pull(.data$symbol)

  purple_cnv_gene %>%
    dplyr::filter(.data$gene %in% genes$symbol) %>%
    dplyr::mutate(
      chromosome = as.factor(.data$chromosome),
      transcriptID = paste0(.data$transcriptId, ".", .data$transcriptVersion),
      minRegStartEnd = paste0(.data$minRegionStart, "-", .data$minRegionEnd),
      minRegSupportStartEndMethod = paste0(.data$minRegionStartSupport, "-", .data$minRegionEndSupport,
                                           " (", .data$minRegionMethod, ")"),
      germDelReg = paste0(.data$germlineHomDeletionRegions, "/", .data$germlineHetToHomDeletionRegions),
      oncogene = .data$gene %in% oncogenes,
      tsgene = .data$gene %in% tsgenes,
      onco_or_ts = dplyr::case_when(
        .data$oncogene & .data$tsgene ~ "onco+ts",
        .data$oncogene ~ "oncogene",
        .data$tsgene ~ "tsgene",
        TRUE ~ "")) %>%
    dplyr::select(.data$gene, minCN = .data$minCopyNumber, maxCN = .data$maxCopyNumber,
                  chrom = .data$chromosome, .data$start, .data$end,
                  chrBand = .data$chromosomeBand, .data$onco_or_ts, .data$transcriptID, .data$minMinorAllelePloidy,
                  somReg = .data$somaticRegions, .data$germDelReg, minReg = .data$minRegions,
                  .data$minRegStartEnd, .data$minRegSupportStartEndMethod)
}

#' Read PURPLE CNV Somatic File
#'
#' Reads the `purple.cnv.somatic.tsv` file, which contains the copy number
#' profile of all (contiguous) segments of the tumor sample
#' (see [this table](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator#copy-number-file)).
#'
#' @param x Path to `purple.cnv.somatic.tsv` file.
#' @param v PURPLE version (default: 2.39). Used to determine the column names.
#'
#' @return The input file as a tibble.
#'
#' @examples
#' x <- system.file("extdata/purple/v2.39/purple.cnv.somatic.tsv", package = "gpgr")
#' p <- read_purple_cnv_somatic(x, "2.39")
#' p
#'
#' @testexamples
#' expect_equal(colnames(p)[ncol(p)], "majorAllelePloidy")
#' expect_error(read_purple_cnv_somatic(x, "2.38"))
#'
#' @export
read_purple_cnv_somatic <- function(x, v = "2.39") {

  config <- function(v) {
    current <- list(
      nm = c("chromosome", "start", "end", "copyNumber", "bafCount", "observedBAF",
             "baf", "segmentStartSupport", "segmentEndSupport", "method",
             "depthWindowCount", "gcContent", "minStart", "maxStart", "minorAllelePloidy",
             "majorAllelePloidy"),
      type = "ciididdcccdddddd")
    l <- list(
      "2.39" = current,
      "2.45" = list(nm = c(current$nm, "foo", "bar"),
                    type = paste0(current$type, "cc")))
    if (numeric_version(v) >= numeric_version("2.45")) {
      return(l[["2.45"]])
    } else if (numeric_version(v) >= numeric_version("2.39")) {
      return(l[["2.39"]])
    } else {
      stop(glue::glue("You've specified PURPLE version {v}, which is older than 2.39.",
                      "Please update to newer PURPLE version."))
    }
  }
  conf <- config(v)
  purple_cnv_somatic <- readr::read_tsv(x, col_types = conf$type)
  assertthat::assert_that(ncol(purple_cnv_somatic) == length(conf$nm))
  assertthat::assert_that(all(colnames(purple_cnv_somatic) == conf$nm))
  purple_cnv_somatic
}

#' Process PURPLE CNV Somatic File for UMCCRISE
#'
#' Processes the `purple.cnv.somatic.tsv` file.
#' and selects columns of interest.
#'
#' @param x Path to `purple.cnv.somatic.tsv` file.
#' @param v PURPLE version (default: 2.39). Used to determine the column names.
#'
#' @return Tibble with more condensed columns.
#'
#' @examples
#' x <- system.file("extdata/purple/v2.39/purple.cnv.somatic.tsv", package = "gpgr")
#' pp <- process_purple_cnv_somatic(x)
#'
#' @testexamples
#' expect_equal(colnames(pp)[ncol(pp)], "GC (windowCount)")
#' expect_error(process_purple_cnv_somatic(x, "2.38"))
#'
#' @export
process_purple_cnv_somatic <- function(x, v = "2.39") {

  purple_cnv_somatic <- read_purple_cnv_somatic(x, v)
  purple_cnv_somatic %>%
    dplyr::mutate(
      Chr = as.factor(.data$chromosome),
      minorAllelePloidy = round(.data$minorAllelePloidy, 1),
      majorAllelePloidy = round(.data$majorAllelePloidy, 1),
      `Ploidy Min+Maj` = paste0(.data$minorAllelePloidy, "+", .data$majorAllelePloidy),
      copyNumber = round(.data$copyNumber, 1),
      bafAdj = round(.data$baf, 2),
      gcContent = round(.data$gcContent, 2),
      `Start/End SegSupport` = paste0(.data$segmentStartSupport, "-", .data$segmentEndSupport),
      `BAF (count)` = paste0(.data$bafAdj, " (", .data$bafCount, ")"),
      `GC (windowCount)` = paste0(.data$gcContent, " (", .data$depthWindowCount, ")")) %>%
    dplyr::select(
      .data$Chr, Start = .data$start, End = .data$end, CN = .data$copyNumber,
      .data$`Ploidy Min+Maj`, .data$`Start/End SegSupport`, Method = .data$method,
      .data$`BAF (count)`, .data$`GC (windowCount)`)
}
