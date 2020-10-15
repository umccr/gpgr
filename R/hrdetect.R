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
#' x <- system.file("extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
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
  assertthat::assert_that(file.exists(x))
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
#' @param genome Human genome version (default: hg38. hg19 means GRCh37).
#'
#' @return Tibble with following BEDPE-like columns:
#' - chrom1, start1, end1
#' - chrom2, start2, end2
#' - sample
#' - strand1, strand2
#'
#'
#' @examples
#' x <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' sv_bedpe <- hrdetect_read_sv_vcf(x, nm = "SAMPLE")
#' head(sv_bedpe)
#'
#' @testexamples
#' expect_equal(nrow(sv_bedpe), 190)
#' expect_equal(colnames(sv_bedpe), c("chrom1", "start1", "end1", "chrom2",
#'              "start2", "end2", "sample", "strand1", "strand2"))
#' @export
hrdetect_read_sv_vcf <- function(x, nm = NULL, genome = "hg38") {

  assertthat::assert_that(file.exists(x))
  assertthat::assert_that(!is.null(nm))
  assertthat::assert_that(genome %in% c("hg19", "hg38", "GRCh37"))
  if (genome == "GRCh37") {
    genome <- "hg19"
  }

  vcf <- VariantAnnotation::readVcf(x, genome)
  gr <- StructuralVariantAnnotation::breakpointRanges(vcf)
  bedpe <- StructuralVariantAnnotation::breakpointgr2bedpe(gr) %>%
    dplyr::mutate_if(is.factor, as.character) %>%
    dplyr::mutate(sample = nm) %>%
    dplyr::select(.data$chrom1, .data$start1, .data$end1,
                  .data$chrom2, .data$start2, .data$end2,
                  .data$sample, .data$strand1, .data$strand2)

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
#' x <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' (cnv <- hrdetect_read_purple_cnv(x))
#'
#' @testexamples
#' expect_equal(colnames(cnv), c("Chromosome", "chromStart", "chromEnd",
#'                               "total.copy.number.inTumour",
#'                               "minor.copy.number.inTumour"))
#'
#' @export
hrdetect_read_purple_cnv <- function(x) {

  assertthat::assert_that(file.exists(x))
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
#' Prepares VCF with SNVs/INDELs for use with HRDetect.
#'
#' @param x Path to VCF with SNVs and INDELs.
#' @param nm Sample name.
#' @param outdir Directory to output analysis results.
#' @param genome Human genome version (default: hg38. hg19 means GRCh37).
#' @param sigsToUse COSMIC signatures to use.
#'
#' @return List with two elements:
#' - snv_results: tibble with exposure score and p-value for chosen signatures.
#' - indel_results: tibble with a summary of the count of indels and their proportion.
#'
#' @examples
#' x <- system.file("extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
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

  assertthat::assert_that(file.exists(x), !is.null(nm), !is.null(outdir),
                          all(sigsToUse %in% 1:30), all(c(3, 8) %in% sigsToUse))
  assertthat::assert_that(genome %in% c("hg19", "hg38", "GRCh37"))
  if (genome == "GRCh37") {
    genome <- "hg19"
  }

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
#' Prepares VCF with SVs for use with HRDetect.
#'
#' @param x Path to VCF with SVs.
#' @param nm Sample name.
#' @param genome Human genome version (default: hg38. hg19 means GRCh37).
#'
#' @return Single-column data.frame (with rownames) with counts for each SV category.
#'
#' @examples
#' x <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' nm <- "SampleA"
#' (d <- hrdetect_prep_sv(x, nm))
#'
#' @testexamples
#' expect_equal(colnames(d), nm)
#' expect_true(inherits(d, "data.frame"))
#'
#' @export
hrdetect_prep_sv <- function(x, nm = NULL, genome = "hg38") {

  assertthat::assert_that(file.exists(x))
  if (genome == "GRCh37") {
    genome <- "hg19"
  }
  sv_bedpe <- hrdetect_read_sv_vcf(x, nm = nm, genome = genome)

  res <- signature.tools.lib::bedpeToRearrCatalogue(sv_bedpe = sv_bedpe)[["rearr_catalogue"]]
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
#' x <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' (l <- hrdetect_prep_cnv(x, nm = "SampleA"))
#'
#' @testexamples
#' expect_equal(colnames(l), c("name", "hrdloh_index"))
#' expect_equal(nrow(l), 1)
#'
#' @export
hrdetect_prep_cnv <- function(x, nm = NULL) {

  assertthat::assert_that(file.exists(x))
  assertthat::assert_that(!is.null(nm))
  cnv <- hrdetect_read_purple_cnv(x)
  cnv_hrd <- signature.tools.lib::ascatToHRDLOH(ascat.data = cnv, SAMPLE.ID = nm)

  tibble::tibble(name = names(cnv_hrd),
                 hrdloh_index = unname(cnv_hrd))
}

#' Run HRDetect via signature.tools.lib
#'
#' Runs HRDetect as described in the
#' [signature.tools.lib repository](https://github.com/Nik-Zainal-Group/signature.tools.lib).
#'
#' @param snvindel_vcf Path to VCF with SNVs and INDELs.
#' @param sv_vcf Path to VCF with SVs.
#' @param cnv_tsv Path to `purple.cnv.somatic.tsv` file.
#' @param nm Sample name.
#' @param genome Human genome version (default: hg38. hg19 means GRCh37).
#' @param snvoutdir Directory to output SNV signature analysis results.
#' @param sigsToUse COSMIC SNV signatures to use.
#' @return Tibble with sample name and HRD probability in first two columns.
#'
#' @examples
#' snvindel_vcf <- system.file(
#'                   "extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz",
#'                   package = "gpgr")
#' sv_vcf <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' cnv_tsv <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' nm <- "SampleA"
#' genome <- "hg38"
#' snvoutdir <- tempdir()
#' (res <- hrdetect_run(nm, snvindel_vcf, sv_vcf, cnv_tsv, genome, snvoutdir))
#'
#' @testexamples
#' expect_equal(colnames(res), c("sample", "Probability", "intercept", "del.mh.prop", "SNV3",
#'                               "SV3", "SV5", "hrd", "SNV8"))
#' expect_true(inherits(res, "data.frame"))
#'
#' @export
hrdetect_run <- function(nm, snvindel_vcf, sv_vcf, cnv_tsv, genome, snvoutdir,
                         sigsToUse = c(1, 2, 3, 5, 6, 8, 13, 17, 18, 20, 26, 30)) {

  assertthat::assert_that(file.exists(snvindel_vcf, sv_vcf, cnv_tsv))
  if (genome == "GRCh37") {
    genome <- "hg19"
  }

  snvindel <- hrdetect_prep_snvindel(snvindel_vcf, nm, genome, snvoutdir, sigsToUse = sigsToUse)
  snv <- snvindel$snv_results %>%
    dplyr::filter(.data$sig %in% c("Signature.3", "Signature.8")) %>%
    tidyr::pivot_wider(id_cols = c("sig", "exposure"),
                       names_from = "sig", values_from = "exposure")

  indel <- snvindel$indel_results$del.mh.prop
  sv <- hrdetect_prep_sv(sv_vcf, nm, genome)
  cnv <- hrdetect_prep_cnv(cnv_tsv, nm)

  tib <- tibble::tibble(
    "del.mh.prop" = indel,
    "SNV3" = snv$Signature.3,
    "SV3" = NA,
    "SV5" = NA,
    "hrd" = cnv$hrdloh_index,
    "SNV8" = snv$Signature.8)
  mat <- as.matrix(tib)
  rownames(mat) <- nm
  res <- signature.tools.lib::HRDetect_pipeline(mat,
                                                genome.v = genome,
                                                SV_catalogues = sv,
                                                nparallel = 2)

  if ("hrdetect_output" %in% names(res)) {
    res <- res[["hrdetect_output"]]
  } else { # no result
    intercept <- matrix(c(NA), dimnames = list(nm, "intercept"))
    Probability <- matrix(c(NA), dimnames = list(nm, "Probability"))
    res <- cbind(intercept, mat, Probability)
  }

  res %>%
    tibble::as_tibble(rownames = "sample", .name_repair = "check_unique") %>%
    dplyr::relocate(.data$Probability, .after = .data$sample)
}
