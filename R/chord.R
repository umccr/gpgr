#' Run CHORD
#'
#' Runs CHORD given SNV and SV VCF files. **NOTE**: make sure you have
#' the BSgenome.Hsapiens.UCSC.hgXX installed.
#'
#' @param vcf.snv Path to VCF containing SNVs and INDELs.
#' @param vcf.sv Path to VCF containing SVs.
#' @param df.sv A data.frame object containing the columns
#'        'SVTYPE' and 'SVLEN' from a Manta SV VCF.
#' @param sample.name Name of sample to use.
#' @param ref.genome Human genome assembly. One of 'hg38' (default), 'hg19' or 'GRCh37'.
#' @param sv.caller manta (default) or gridss.
#' @param ... Other arguments to be passed to [CHORD::extractSigsChord()].
#'
#' @examples
#'
#' snv <- system.file("extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' sv <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' chord_res <- chord_run(vcf.snv = snv, df.sv = gpgr:::chord_mantavcf2df(sv), sample.name = "foo")
#' # chord_res <- chord_run(vcf.snv = snv, vcf.sv = sv, sample.name = "foo") # a bit slower
#'
#' @testexamples
#'
#' expect_equal(length(chord_res), 2)
#' expect_equal(names(chord_res), c("contexts", "prediction"))
#' expect_equal(chord_res[["prediction"]][, 1, drop = TRUE], "foo")
#'
#'
#' @return List with extracted signatures and HRD prediction.
#'
#' @export
chord_run <- function(vcf.snv = NULL, vcf.sv = NULL, df.sv = NULL, sample.name = NULL, ref.genome = "hg38", sv.caller = "manta", ...) {
  avail_genomes <- c("hg19", "hg38", "GRCh37")
  g <- NULL
  if (ref.genome == "hg38") {
    g <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (ref.genome == "hg19") {
    g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else if (ref.genome == "GRCh37") {
    g <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  } else {
    stop(glue::glue("ref.genome should be one of {paste(avail_genomes, collapse = ', ')}"))
  }

  contexts <- CHORD::extractSigsChord(
    vcf.snv = vcf.snv,
    vcf.sv = vcf.sv,
    df.sv = df.sv,
    sv.caller = sv.caller,
    sample.name = sample.name,
    ref.genome = g,
    ...
  )

  prediction <- CHORD::chordPredict(features = contexts,
                                    rf.model = CHORD::CHORD,
                                    do.bootstrap = TRUE, verbose = FALSE)

  list(
    contexts = contexts,
    prediction = prediction
  )
}

#' Convert Manta VCF to data.frame
#'
#' Converts a Manta VCF to a data.frame for processing with CHORD.
#'
#' @param in_vcf Manta VCF.
#'
#' @return Tibble with two columns: sv_type and sv_len (INFO/SVTYPE and INFO/SVLEN from VCF).
#'
#' @examples
#' in_vcf <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' d <- chord_mantavcf2df(in_vcf)
#'
#' @testexamples
#' expect_equal(d$sv_len[1], "-108")
#' expect_equal(d$sv_type[1], "DEL")
#'
#' @export
chord_mantavcf2df <- function(in_vcf) {
  assertthat::assert_that(grepl("vcf$", in_vcf) | grepl("vcf\\.gz$", in_vcf))
  d <- bedr::read.vcf(in_vcf, split.info = TRUE, verbose = FALSE)
  tibble::tibble(sv_type = d$vcf$SVTYPE,
                 sv_len = d$vcf$SVLEN)
}
