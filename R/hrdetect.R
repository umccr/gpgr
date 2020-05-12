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
#' (hrdetect_read_snvindel_vcf(x))
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
#' @param hg_version Human genome version (default: hg38).
#'
#' @return Tibble with following BEDPE-like columns:
#' - chrom1
#' - start1
#' - end1
#' - chrom2
#' - start2
#' - end2
#' - sample
#' - strand1
#' - strand2
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
hrdetect_read_sv_vcf <- function(x, nm = NULL, hg_version = "hg38") {

  assertthat::assert_that(!is.null(nm))

  vcf <- VariantAnnotation::readVcf(x, hg_version)
  gr <- StructuralVariantAnnotation::breakpointRanges(vcf)
  bedpe <- StructuralVariantAnnotation::breakpointgr2bedpe(gr) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate(sample = nm) %>%
    dplyr::select(.data$chrom1:.data$end2, .data$sample, .data$strand1, .data$strand2)

  bedpe

}
