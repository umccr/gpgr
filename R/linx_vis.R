#' Read LINX VisCopyNumber File
#'
#' Reads the `linx.vis_copy_number.tsv` file.
#'
#' @param x Path to `linx.vis_copy_number.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_copy_number.tsv.gz", package = "gpgr")
#' (l <- linx_viscopynumber_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "BAF")
#'
#' @export
linx_viscopynumber_read <- function(x) {
  nm <- c(
    SampleId = "c", Chromosome = "c", Start = "i", End = "i",
    CopyNumber = "d", BAF = "d"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d |>
    dplyr::mutate(Chrom = sub("chr", "", .data$Chromosome)) |>
    dplyr::select(.data$Chrom, .data$Start, .data$End, CN = .data$CopyNumber, .data$BAF) |>
    dplyr::arrange(mixedrank(.data$Chrom)) |>
    dplyr::mutate(Chrom = as.factor(.data$Chrom))
}

#' Read LINX VisFusion File
#'
#' Reads the `linx.vis_fusion.tsv` file.
#'
#' @param x Path to `linx.vis_fusion.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_fusion.tsv.gz", package = "gpgr")
#' (l <- linx_visfusion_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "FusedExonDown")
#'
#' @export
linx_visfusion_read <- function(x) {
  nm <- c(
    SampleId = "c", ClusterId = "d", Reportable = "c",
    GeneNameUp = "c", TranscriptUp = "c", ChrUp = "c",
    PosUp = "d", StrandUp = "d", RegionTypeUp = "c",
    FusedExonUp = "d", GeneNameDown = "c", TranscriptDown = "c",
    ChrDown = "c", PosDown = "d", StrandDown = "d",
    RegionTypeDown = "c", FusedExonDown = "d"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d
}

#' Read LINX VisGeneExon File
#'
#' Reads the `linx.vis_gene_exon.tsv` file.
#'
#' @param x Path to `linx.vis_gene_exon.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_gene_exon.tsv.gz", package = "gpgr")
#' (l <- linx_visgeneexon_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "ExonEnd")
#'
#' @export
linx_visgeneexon_read <- function(x) {
  nm <- c(
    SampleId = "c", ClusterId = "d", Gene = "c", Transcript = "c",
    Chromosome = "c", AnnotationType = "c", ExonRank = "d",
    ExonStart = "d", ExonEnd = "d"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d |>
    dplyr::select(-.data$SampleId) |>
    dplyr::mutate(Chromosome = sub("chr", "", .data$Chromosome)) |>
    dplyr::arrange(.data$ClusterId, .data$Gene, .data$Chromosome) |>
    dplyr::mutate(ClusterId = factor(.data$ClusterId)) |>
    dplyr::rename(Chrom = .data$Chromosome)
}

#' Read LINX VisProteinDomain File
#'
#' Reads the `linx.vis_protein_domain.tsv` file.
#'
#' @param x Path to `linx.vis_protein_domain.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_protein_domain.tsv.gz", package = "gpgr")
#' (l <- linx_visproteindomain_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "Info")
#'
#' @export
linx_visproteindomain_read <- function(x) {
  nm <- c(
    SampleId = "c", ClusterId = "d", Transcript = "c", Chromosome = "c",
    Start = "i", End = "i", Info = "c"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d |>
    dplyr::select(-.data$SampleId) |>
    dplyr::mutate(Chromosome = sub("chr", "", .data$Chromosome)) |>
    dplyr::arrange(.data$ClusterId, .data$Chromosome, .data$Start) |>
    dplyr::mutate(ClusterId = factor(.data$ClusterId)) |>
    dplyr::rename(Chrom = .data$Chromosome)
}

#' Read LINX VisSegments File
#'
#' Reads the `linx.vis_segments.tsv` file.
#'
#' @param x Path to `linx.vis_segments.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_segments.tsv.gz", package = "gpgr")
#' (l <- linx_vissegments_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "InDoubleMinute")
#'
#' @export
linx_vissegments_read <- function(x) {
  nm <- c(
    SampleId = "c", ClusterId = "d", ChainId = "d", Chromosome = "c",
    PosStart = "c", PosEnd = "c", LinkPloidy = "d", InDoubleMinute = "c"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d |>
    dplyr::select(-.data$SampleId) |>
    dplyr::mutate(Chromosome = sub("chr", "", .data$Chromosome)) |>
    dplyr::arrange(.data$ClusterId, .data$ChainId, .data$Chromosome) |>
    dplyr::mutate(ClusterId = factor(.data$ClusterId))
}

#' Read LINX VisSvData File
#'
#' Reads the `linx.vis_sv_data.tsv` file.
#'
#' @param x Path to `linx.vis_sv_data.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_sv_data.tsv.gz", package = "gpgr")
#' (l <- linx_vissvdata_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "InDoubleMinute")
#'
#' @export
linx_vissvdata_read <- function(x) {
  nm <- c(
    SampleId = "c", ClusterId = "d", ChainId = "d", SvId = "d",
    Type = "c", ResolvedType = "c", IsSynthetic = "c", ChrStart = "c",
    ChrEnd = "c", PosStart = "d", PosEnd = "d", OrientStart = "d",
    OrientEnd = "d", InfoStart = "c", InfoEnd = "c",
    JunctionCopyNumber = "d", InDoubleMinute = "c"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d |>
    dplyr::select(-.data$SampleId) |>
    dplyr::arrange(.data$ClusterId, .data$ChainId) |>
    dplyr::mutate(ClusterId = factor(.data$ClusterId)) |>
    dplyr::select(
      .data$ClusterId, .data$ChrStart, .data$ChrEnd,
      .data$PosStart, .data$PosEnd, dplyr::everything()
    )
}
