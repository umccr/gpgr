#' Read VCF with SNVs/INDELs for use with HRDetect
#'
#' Reads a VCF with SNVs/INDELs for use with HRDetect.
#'
#' @param x Path to VCF.
#'
#' @return List containing CHROM, POS, REF and ALT columns
#'         for SNVs and INDELs in separate tibbles.
#'
#' @examples
#' x <- system.file("extdata/umccrise/v0.18/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' (l <- hrdetect_read_snvindel_vcf(x))
#'
#' @testexamples
#' expect_equal(length(l), 2)
#' expect_equal(names(l), c("snv", "indel"))
#' expect_equal(colnames(l$snv), c("chr", "position", "REF", "ALT"))
#' expect_equal(colnames(l$indel), c("chr", "position", "REF", "ALT"))
#'
#' @export
hrdetect_read_snvindel_vcf <- function(x) {
  ALLOWED_BASES <- c("A", "C", "G", "T")
  d <- x %>%
    readr::read_tsv(
      comment = "##",
      col_types = readr::cols_only("#CHROM" = "c", "POS" = "i", "REF" = "c", "ALT" = "c")) %>%
    dplyr::mutate(vartype = dplyr::case_when(
      .data$REF %in% ALLOWED_BASES & .data$ALT %in% ALLOWED_BASES ~ "SNV",
      TRUE ~ "INDEL")) %>%
    dplyr::rename(chr = "#CHROM", position = "POS")

  snv <- d %>% dplyr::filter(.data$vartype == "SNV") %>% dplyr::select(-.data$vartype)
  indel <- d %>% dplyr::filter(.data$vartype == "INDEL") %>% dplyr::select(-.data$vartype)

  list(
    snv = snv,
    indel = indel
  )
}

#' Read VCF with SVs for use with HRDetect
#'
#' Reads a VCF with SVs for use with HRDetect.
#'
#' @param x Path to VCF.
#' @param nm Sample name.
#' @param genome Human genome version (default: hg38).
#'
#' @return Tibble with following BEDPE-like columns:
#' - chrom1, start1, end1
#' - chrom2, start2, end2
#' - sample
#' - strand1, strand2
#'
#'
#' @examples
#' x <- system.file("extdata/umccrise/v0.18/sv/manta.vcf.gz", package = "gpgr")
#' sv_bedpe <- hrdetect_read_sv_vcf(x, nm = "SAMPLE")
#' head(sv_bedpe)
#'
#' @testexamples
#' expect_equal(nrow(sv_bedpe), 190)
#' expect_equal(colnames(sv_bedpe), c("chrom1", "start1", "end1", "chrom2",
#'              "start2", "end2", "sample", "strand1", "strand2"))
#' @export
hrdetect_read_sv_vcf <- function(x, nm = NULL, genome = "hg38") {

  assertthat::assert_that(!is.null(nm))

  vcf <- VariantAnnotation::readVcf(x, genome)
  gr <- StructuralVariantAnnotation::breakpointRanges(vcf)
  bedpe <- StructuralVariantAnnotation::breakpointgr2bedpe(gr) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate(sample = nm) %>%
    dplyr::select(.data$chrom1:.data$end2, .data$sample, .data$strand1, .data$strand2)

  bedpe

}

#' Read PURPLE Somatic CNVs for HRDetect
#'
#' Reads PURPLE somatic CNVs for use with HRDetect.
#'
#' @param x Path to `purple.cnv.somatic.tsv` file.
#'
#' @return Tibble containing following columns:
#' - chromosome, start, end
#' - copyNumber (total)
#' - minorAllelePloidy
#'
#' @examples
#' x <- system.file("extdata/purple/v2.39/purple.cnv.somatic.tsv", package = "gpgr")
#' (cnv <- hrdetect_read_purple_cnv(x))
#'
#' @testexamples
#' expect_equal(colnames(cnv), c("Chromosome", "chromStart", "chromEnd",
#'                               "total.copy.number.inTumour",
#'                               "minor.copy.number.inTumour"))
#'
#' @export
hrdetect_read_purple_cnv <- function(x) {

  cnv <- readr::read_tsv(x,
                         col_types = readr::cols_only(
                           "chromosome" = "c", "start" = "i", "end" = "i",
                           "copyNumber" = "d", "minorAllelePloidy" = "d"))

  cnv %>%
    dplyr::rename(
      Chromosome = .data$chromosome,
      chromStart = .data$start,
      chromEnd = .data$end,
      total.copy.number.inTumour = .data$copyNumber,
      minor.copy.number.inTumour = .data$minorAllelePloidy,
    ) %>%
    dplyr::mutate(Chromosome = sub("chr", "", .data$Chromosome))
}

#' Prepare VCF with SNVs/INDELs for use with HRDetect
#'
#' @param x Path to VCF with SNVs and INDELs.
#' @param nm Sample name.
#' @param genome Genome assembly version.
#' @param outdir Directory to output analysis results.
#' @param sigsToUse COSMIC signatures to use.
#'
#' @return List with two elements:
#' - snv_results: tibble with exposure score and p-value for chosen signatures.
#' - indel_results: tibble with a summary of the count of indels and their proportion.
#'
#' @examples
#' x <- system.file("extdata/umccrise/v0.18/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' (l <- hrdetect_prep_snvindel(x, nm = "sampleA", outdir = tempdir()))
#'
#' @testexamples
#' expect_equal(c("snv_results", "indel_results"), names(l))
#' expect_equal(c("sig", "exposure", "pvalue"), colnames(l[["snv_results"]]))
#' expect_equal(colnames(l[["indel_results"]])[c(1, 7)], c("sample", "del.mh.prop"))
#'
#' @export
hrdetect_prep_snvindel <- function(x, nm = NULL, genome = "hg38", outdir = NULL,
                                   sigsToUse = c(1, 2, 3, 5, 6, 8, 13, 17, 18, 20, 26, 30)) {

  assertthat::assert_that(!is.null(nm), !is.null(outdir), genome %in% c("hg19", "hg38"),
                          all(sigsToUse %in% 1:30))

  # must end in /
  outdir <- ifelse(!grepl("/$", outdir), paste0(outdir, "/"), outdir)

  snvindel_tabs <- hrdetect_read_snvindel_vcf(x)

  ##--- SNVs ---##
  snv_catalogue <- signature.tools.lib::tabToSNVcatalogue(
    subs = snvindel_tabs[["snv"]],
    genome.v = genome)[["catalogue"]]

  subs_fit_res <- signature.tools.lib::SignatureFit_withBootstrap_Analysis(
    outdir = outdir,
    cat = snv_catalogue,
    signature_data_matrix = signature.tools.lib::COSMIC30_subs_signatures[, sigsToUse],
    type_of_mutations = "subs",
    nboot = 100,
    nparallel = 2)

  snv_exp <- subs_fit_res$E_median_filtered %>%
    tibble::as_tibble(rownames = "sig") %>%
    dplyr::rename(exposure = .data$catalogue)

  snv_pval <- subs_fit_res$E_p.values %>%
    tibble::as_tibble(rownames = "sig") %>%
    dplyr::rename(pvalue = .data$catalogue)

  snv_results <- snv_exp %>%
    dplyr::left_join(snv_pval, by = "sig")

  ##--- INDELs ---##
  indel_count_proportion <- signature.tools.lib::tabToIndelsClassification(
    indel.data = snvindel_tabs[["indel"]],
    sampleID = nm,
    genome.v = genome)[["count_proportion"]] %>%
    tibble::as_tibble()

  list(
    snv_results = snv_results,
    indel_results = indel_count_proportion
  )
}

#' Prepare VCF with SVs for use with HRDetect
#'
#' @param x Path to VCF with SVs.
#' @param nm Sample name.
#' @param genome Human genome version (default: hg38).
#'
#' @return Tibble with counts for each SV category.
#'
#' @examples
#' x <- system.file("extdata/umccrise/v0.18/sv/manta.vcf.gz", package = "gpgr")
#' (l <- hrdetect_prep_sv(x, nm = "SampleA"))
#'
#' @testexamples
#' expect_equal(colnames(l), c("sv_category", "count"))
#'
#' @export
hrdetect_prep_sv <- function(x, nm = NULL, genome = "hg38") {
  sv_bedpe <- hrdetect_read_sv_vcf(x, nm = nm, genome = genome)
  res <- signature.tools.lib::bedpeToRearrCatalogue(sv_bedpe = sv_bedpe)[["rearr_catalogue"]] %>%
    tibble::as_tibble(rownames = "sv_category") %>%
    dplyr::rename(count = nm)
  res
}

#' Prepare PURPLE Somatic CNVs for HRDetect
#'
#' Prepares PURPLE somatic CNVs for use with HRDetect.
#'
#' @param x Path to `purple.cnv.somatic.tsv` file.
#' @param nm Sample name.
#'
#' @return Tibble with sample name and HRD-LOH index.
#'
#' @examples
#' x <- system.file("extdata/purple/v2.39/purple.cnv.somatic.tsv", package = "gpgr")
#' (l <- hrdetect_prep_cnv(x, nm = "SampleA"))
#'
#' @testexamples
#' expect_equal(colnames(l), c("name", "hrdloh_index"))
#' expect_equal(nrow(l), 1)
#'
#' @export
hrdetect_prep_cnv <- function(x, nm = NULL) {

  assertthat::assert_that(!is.null(nm))
  cnv <- hrdetect_read_purple_cnv(x)
  cnv_hrd <- signature.tools.lib::ascatToHRDLOH(ascat.data = cnv, SAMPLE.ID = nm)

  tibble::tibble(name = names(cnv_hrd),
                 hrdloh_index = unname(cnv_hrd))
}

