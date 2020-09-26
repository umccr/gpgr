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
#' @param ... Other arguments to be passed to `CHORD::extractSigsChord()`.
#'
#' @examples
#' \dontrun{
#' snv <- system.file("extdata/umccrise/v0.18/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' sv <- system.file("extdata/umccrise/v0.18/sv/manta.vcf.gz", package = "gpgr")
#' chord_res <- run_chord(vcf.snv = snv, vcf.sv = sv, sample.name = "foo")
#' chord_res2 <- run_chord(vcf.snv = snv, df.sv = gpgr:::chord_mantavcf2df(sv), sample.name = "foo")
#' }
#'
#'
#' @return List with extracted signatures and HRD prediction.
#'
#' @export
run_chord <- function(vcf.snv = NULL, vcf.sv = NULL, df.sv = NULL, sample.name = NULL, ref.genome = "hg38", sv.caller = "manta", ...) {
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

chord_mantavcf2df <- function(in_vcf) {
  d <- bedr::read.vcf(in_vcf, split.info = TRUE, verbose = FALSE)
  tibble::tibble(sv_type = d$vcf$SVTYPE,
                 sv_len = d$vcf$SVLEN)
}
