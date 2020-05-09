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

