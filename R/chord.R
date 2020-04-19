#' Run CHORD
#'
#' Runs CHORD given SNV and SV VCF files.
#'
#' @param snv Path to SNV VCF.
#' @param sv Path to SV VCF.
#' @param sample Name of sample to use.
#' @param genome Human genome assembly. One of 'hg38' (default) or 'hg19'.
#'
#' @return List with extracted signatures and HRD prediction.
#' @export
#' @import BSgenome.Hsapiens.UCSC.hg38, mutSigExtractor, CHORD
#'
#' @examples
run_chord <- function(snv, sv, sample, genome = "hg38") {
  assertthat::assert_that(all(file.exists(c(snv, sv))))
  assertthat::assert_that(genome %in% c("hg19", "hg38"))
  ref_genome <- paste0("BSgenome.Hsapiens.UCSC.", genome)

  contexts <- CHORD::extractSigsChord(
    vcf.snv = snv,
    vcf.sv = sv,
    sv.caller='manta',
    sample.name = sample,
    ref.genome = ref_genome
  )

  prediction <- CHORD::chordPredict(contexts, do.bootstrap = TRUE, verbose = FALSE)

  list(
    contexts = contexts,
    prediction = prediction
  )

}
