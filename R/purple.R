#' Read PURPLE CNV Gene File
#'
#' Reads the `purple.cnv.gene.tsv` file, which summarises copy number
#' alterations of each gene in the HMF panel
#' (see [this table](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator#gene-copy-number-file)).
#'
#' @param x Path to `purple.cnv.gene.tsv` file.
#'
#' @return The input file as a tibble.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.gene.tsv", package = "gpgr")
#' (p <- read_purple_cnv_gene(x))
#'
#' @testexamples
#' expect_equal(colnames(p)[ncol(p)], "minMinorAllelePloidy")
#'
#' @export
read_purple_cnv_gene <- function(x) {

  nm <- c("chromosome" = "c", "start" = "i", "end" = "i", "gene" = "c",
          "minCopyNumber" = "d", "maxCopyNumber" = "d",
          "unused" = "c", "somaticRegions" = "d", "germlineHomDeletionRegions" = "d",
          "germlineHetToHomDeletionRegions" = "d",
          "transcriptId" = "c", "transcriptVersion" = "c", "chromosomeBand" = "c",
          "minRegions" = "d", "minRegionStart" = "i", "minRegionEnd" = "i",
          "minRegionStartSupport" = "c", "minRegionEndSupport" = "c",
          "minRegionMethod" = "c", "minMinorAllelePloidy" = "d")

  ctypes <- paste(nm, collapse = "")
  purple_cnv_gene <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(purple_cnv_gene) == length(nm))
  assertthat::assert_that(all(colnames(purple_cnv_gene) == names(nm)))
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
#'
#' @return List with two elements:
#' * `tab`: Tibble filtered to genes found in  `g`.
#' * `descr`: Description of tibble columns.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.gene.tsv", package = "gpgr")
#' g <- system.file("extdata/ref/umccr_cancer_genes_2019-03-20.tsv", package = "gpgr")
#' (pp <- process_purple_cnv_gene(x, g))
#'
#' @testexamples
#' expect_equal(colnames(pp$tab)[ncol(pp$tab)], "minRegSupportStartEndMethod")
#'
#' @export
process_purple_cnv_gene <- function(x, g = NULL) {
  purple_cnv_gene <- read_purple_cnv_gene(x)
  if (is.null(g)) {
    g <- system.file("extdata/ref/umccr_cancer_genes_2019-03-20.tsv", package = "gpgr")
  }
  genes <- readr::read_tsv(g, col_types = readr::cols(symbol = "c", oncogene = "l", tumorsuppressor = "l")) %>%
    dplyr::select(.data$symbol, .data$oncogene, .data$tumorsuppressor)
  oncogenes <- genes %>% dplyr::filter(.data$oncogene) %>% dplyr::pull(.data$symbol)
  tsgenes <- genes %>% dplyr::filter(.data$tumorsuppressor) %>% dplyr::pull(.data$symbol)

  purple_cnv_gene <- purple_cnv_gene %>%
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
                  chrBand = .data$chromosomeBand, .data$onco_or_ts,
                  .data$transcriptID, .data$minMinorAllelePloidy,
                  somReg = .data$somaticRegions, .data$germDelReg, minReg = .data$minRegions,
                  .data$minRegStartEnd, .data$minRegSupportStartEndMethod)

  descr <- dplyr::tribble(
    ~Column, ~Description,
    "gene", "Name of gene",
    "minCN/maxCN", "Min/Max copy number found in gene exons",
    "chrom/start/end", "Chromosome/start/end location of gene transcript",
    "chrBand", "Chromosome band of the gene",
    "onco_or_ts", "oncogene ('oncogene'), tumor suppressor ('tsgene'), or both ('onco+ts'), as reported by [Cancermine](https://github.com/jakelever/cancermine)",
    "transcriptID", "Ensembl transcript ID (dot version)",
    "minMinorAllelePloidy", "Minimum allele ploidy found over the gene exons - useful for identifying LOH events",
    "somReg (somaticRegions)", "Count of somatic copy number regions this gene spans",
    "germDelReg (germlineHomDeletionRegions / germlineHetToHomDeletionRegions)", "Number of regions spanned by this gene that are (homozygously deleted in the germline / both heterozygously deleted in the germline and homozygously deleted in the tumor)",
    "minReg (minRegions)", "Number of somatic regions inside the gene that share the min copy number",
    "minRegStartEnd", "Start/End base of the copy number region overlapping the gene with the minimum copy number",
    "minRegSupportStartEndMethod", "Start/end support of the CN region overlapping the gene with the min CN (plus determination method)")

  list(tab = purple_cnv_gene,
       descr = descr)
}

#' Read PURPLE CNV Somatic File
#'
#' Reads the `purple.cnv.somatic.tsv` file, which contains the copy number
#' profile of all (contiguous) segments of the tumor sample
#' (see [this table](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator#copy-number-file)).
#'
#' @param x Path to `purple.cnv.somatic.tsv` file.
#'
#' @return The input file as a tibble.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' (p <- read_purple_cnv_somatic(x))
#'
#' @testexamples
#' expect_equal(colnames(p)[ncol(p)], "majorAllelePloidy")
#'
#' @export
read_purple_cnv_somatic <- function(x) {
  nm <- c("chromosome" = "c", "start" = "i", "end" = "i",
          "copyNumber" = "d", "bafCount" = "d", "observedBAF" = "d",
          "baf" = "d", "segmentStartSupport" = "c", "segmentEndSupport" = "c",
          "method" = "c", "depthWindowCount" = "i", "gcContent" = "d",
          "minStart" = "i", "maxStart" = "i", "minorAllelePloidy" = "d",
          "majorAllelePloidy" = "d")
  ctypes <- paste(nm, collapse = "")
  purple_cnv_somatic <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(purple_cnv_somatic) == length(nm))
  assertthat::assert_that(all(colnames(purple_cnv_somatic) == names(nm)))
  purple_cnv_somatic
}

#' Process PURPLE CNV Somatic File for UMCCRISE
#'
#' Processes the `purple.cnv.somatic.tsv` file.
#' and selects columns of interest.
#'
#' @param x Path to `purple.cnv.somatic.tsv` file.
#'
#' @return List with two elements:
#' * `tab`: Tibble with more condensed columns.
#' * `descr`: Description of tibble columns.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' (pp <- process_purple_cnv_somatic(x))
#'
#' @testexamples
#' expect_equal(colnames(pp$tab)[ncol(pp$tab)], "GC (windowCount)")
#'
#' @export
process_purple_cnv_somatic <- function(x) {

  purple_cnv_somatic <- read_purple_cnv_somatic(x)
  purple_cnv_somatic <- purple_cnv_somatic %>%
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


  descr <- dplyr::tribble(
    ~Column, ~Description,
    "Chr/Start/End", "Coordinates of copy number segment",
    "CN", "Fitted absolute copy number of segment adjusted for purity and ploidy",
    "Ploidy Min+Maj", "Ploidy of minor + major allele adjusted for purity",
    "Start/End SegSupport", paste0("Type of SV support for the CN breakpoint at ",
                                   "start/end of region. Allowed values: ",
                                   "CENTROMERE, TELOMERE, INV, DEL, DUP, BND (translocation), ",
                                   "SGL (single breakend SV support), NONE (no SV support for CN breakpoint), ",
                                   "MULT (multiple SV support at exact breakpoint)"),
    "Method", paste0("Method used to determine the CN of the region. Allowed values: ",
                     "BAF_WEIGHTED (avg of all depth windows for the region), ",
                     "STRUCTURAL_VARIANT (inferred using ploidy of flanking SVs), ",
                     "LONG_ARM (inferred from the long arm), GERMLINE_AMPLIFICATION ",
                     "(inferred using special logic to handle regions of germline amplification)"),
    "BAF (count)", "Tumor BAF after adjusted for purity and ploidy (Count of AMBER baf points covered by this segment)",
    "GC (windowCount)", "Proportion of segment that is G or C (Count of COBALT windows covered by this segment)"
  )

  list(tab = purple_cnv_somatic,
       descr = descr)
}

#' Read PURPLE CNV Germline File
#'
#' Reads the `purple.cnv.germline.tsv` file.
#'
#' @param x Path to `purple.cnv.germline.tsv` file.
#'
#' @return The input file as a tibble.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.germline.tsv", package = "gpgr")
#' (p <- read_purple_cnv_germline(x))
#'
#' @testexamples
#' expect_equal(colnames(p)[ncol(p)], "majorAllelePloidy")
#'
#' @export
read_purple_cnv_germline <- function(x) {
  # as of PURPLE v2.39, germline and somatic files have same columns.
  purple_cnv_germline <- read_purple_cnv_somatic(x)
  purple_cnv_germline
}

#' Process PURPLE CNV germline File for UMCCRISE
#'
#' Processes the `purple.cnv.germline.tsv` file.
#' and selects columns of interest.
#'
#' @param x Path to `purple.cnv.germline.tsv` file.
#'
#' @return List with two elements:
#' * `tab`: Tibble with more condensed columns.
#' * `descr`: Description of tibble columns.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.germline.tsv", package = "gpgr")
#' (pp <- process_purple_cnv_germline(x))
#'
#' @testexamples
#' expect_equal(colnames(pp$tab)[ncol(pp$tab)], "GC (windowCount)")
#'
#' @export
process_purple_cnv_germline <- function(x) {
  # as of PURPLE v2.39, germline and somatic files have same columns.
  processed_purple_cnv_germline <- process_purple_cnv_somatic(x)
  processed_purple_cnv_germline
}

#' Read PURPLE version file
#'
#' Reads the `purple.version` file containing the PURPLE version and build date.
#'
#' @param x Path to the `purple.version` file.
#'
#' @return A list with the elements:
#' * `version`: The version of PURPLE used.
#' * `build_date`: The build date of the PURPLE version.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.version", package = "gpgr")
#' (v <- read_purple_version(x))
#'
#' @testexamples
#' expect_equal(length(v), 2)
#' expect_equal(names(v), c("version", "build_date"))
#' expect_equal(v$version, "2.39")
#'
#' @export
read_purple_version <- function(x) {
  tab <- readr::read_delim(x, delim = "=", col_names = c("key", "value"), col_types = "cc")
  assertthat::assert_that(nrow(tab) == 2, all(tab$key == c("version", "build.date")))

  list(version = tab$value[tab$key == "version"],
       build_date = tab$value[tab$key == "build.date"])
}

#' Read PURPLE QC file
#'
#' Reads the `purple.qc` file.
#'
#' @param x Path to the `purple.qc` file.
#'
#' @return A named character vector containing QC values for several PURPLE
#'         metrics.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.qc", package = "gpgr")
#' (p <- purple_qc_read(x))
#'
#' @testexamples
#' expect_true(p[1] == "PASS")
#'
#' @export
purple_qc_read <- function(x) {
  purple_qc <-
    readr::read_tsv(x, col_names = c("key", "value"), col_types = "cc") %>%
    dplyr::mutate(value = toupper(.data$value))

  nm <- c("QCStatus", "SegmentPass", "GenderPass", "DeletedGenesPass",
          "SegmentScore", "UnsupportedSegments", "Ploidy", "AmberGender",
          "CobaltGender", "DeletedGenes")

  assertthat::assert_that(all(purple_qc$key == nm))
  structure(purple_qc$value, names = purple_qc$key)
}

#' Read PURPLE Purity file
#'
#' Reads the `purple.purity.tsv` file containing a summary of the purity fit.
#'
#' @param x Path to the `purple.purity.tsv` file.
#'
#' @return The input file as a tibble with a description of each metric.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.purity.tsv", package = "gpgr")
#' (p <- purple_purity_read(x))
#'
#' @testexamples
#' expect_equal(p[1, "Column", drop = TRUE], "purity")
#' expect_equal(p[nrow(p), "Column", drop = TRUE], "tmbStatus")
#'
#' @export
purple_purity_read <- function(x) {
  tab <- dplyr::tribble(
    ~Column, ~Description, ~Type,
    "purity", "Purity of tumor in the sample.", "d",
    "normFactor", "Factor to convert tumor ratio to copy number. Lower number implies higher ploidy.", "d",
    "score", "Score of fit. Lower is better.", "d",
    "diploidProportion", "Proportion of CN regions that have 1 (+- 0.2) minor and major allele.", "d",
    "ploidy", "Average ploidy of tumor after adjusting for purity.", "d",
    "gender", "One of MALE, FEMALE or MALE_KLINEFELTER.", "c",
    "status", "One of NORMAL, HIGHLY_DIPLOID, SOMATIC or NO_TUMOR.", "c",
    "polyclonalProportion", "Proportion of CN regions that are more than 0.25 from a whole CN.", "d",
    "minPurity", "Minimum purity with score within 10% of best.", "d",
    "maxPurity", "Maximum purity with score within 10% of best.", "d",
    "minPloidy", "Minimum ploidy with score within 10% of best.", "d",
    "maxPloidy", "Maximum ploidy with score within 10% of best.", "d",
    "minDiploidProportion", "Min diploidProportion.", "d",
    "maxDiploidProportion", "Max diploidProportion.", "d",
    "version", "PURPLE Version.", "c",
    "somaticPenalty", paste("Penalty from somatic variants with implied variant",
                            "copy numbers that are inconsistent with the minor and major allele copy number."), "d",
    "wholeGenomeDuplication", "True if more than 10 autosomes have major allele copy number > 1.5.", "c",
    "msIndelsPerMb", "Microsatellite indels per mega base.", "d",
    "msStatus", "Microsatellite status. One of MSI, MSS or UNKNOWN if somatic variants not supplied.", "c",
    "tml", "Tumor mutational load (# of missense variants in sample).", "d",
    "tmlStatus", "Tumor mutational load status. One of HIGH, LOW or UNKNOWN if somatic variants not supplied.", "c",
    "tmbPerMb", "Tumor mutational burden (#passing variants per Mb) per mega base.", "d",
    "tmbStatus", paste("Tumor mutational burden status. One of HIGH, LOW or UNKNOWN if",
                       "somatic variants not supplied. High = > 10 pass variants per Mb."), "c")

  ctypes <- paste(tab$Type, collapse = "")
  purple_purity <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(purple_purity) == nrow(tab))
  assertthat::assert_that(all(colnames(purple_purity) == tab$Column))

  purple_purity %>%
    dplyr::mutate(
      dplyr::across(tidyselect::vars_select_helpers$where(is.numeric), round, 2),
      dplyr::across(dplyr::everything(), as.character)) %>%
    tidyr::pivot_longer(dplyr::everything(), names_to = "Column", values_to = "Value") %>%
    dplyr::left_join(tab, by = "Column") %>%
    dplyr::select(.data$Column, .data$Value, .data$Description)
}
