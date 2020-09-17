#' Run CHORD
#'
#' Runs CHORD given SNV and SV VCF files. **NOTE**: make sure you have
#' the BSgenome.Hsapiens.UCSC.hgXX installed.
#'
#' @param snv Path to SNV VCF.
#' @param sv Path to SV VCF.
#' @param sample Name of sample to use.
#' @param genome Human genome assembly. One of 'hg38' (default) or 'hg19'.
#'
#' @examples
#' \dontrun{
#' snv <- system.file("extdata/umccrise/v0.18/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' sv <- system.file("extdata/umccrise/v0.18/sv/manta.vcf.gz", package = "gpgr")
#' chord_res <- run_chord(snv = snv, sv = sv, sample = "foo")
#' }
#'
#'
#' @return List with extracted signatures and HRD prediction.
#'
#' @export
run_chord <- function(snv, sv, sample, genome = "hg38") {
  assertthat::assert_that(all(file.exists(c(snv, sv))))
  assertthat::assert_that(genome %in% c("hg19", "hg38", "GRCh37"))
  ref_genome <- NULL
  if (genome == "hg38") {
    ref_genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
  } else if (genome == "hg19") {
    ref_genome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  } else if (genome == "GRCh37") {
    ref_genome <- BSgenome.Hsapiens.1000genomes.hs37d5::BSgenome.Hsapiens.1000genomes.hs37d5
  }

  contexts <- CHORD::extractSigsChord(
    vcf.snv = snv,
    vcf.sv = sv,
    sv.caller='manta',
    sample.name = sample,
    ref.genome = ref_genome
  )

  prediction <- CHORD::chordPredict(features = contexts,
                                    rf.model = CHORD::CHORD,
                                    do.bootstrap = TRUE, verbose = FALSE)

  list(
    contexts = contexts,
    prediction = prediction
  )

}
