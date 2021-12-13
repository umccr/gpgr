#' Read LINX SVS File
#'
#' Reads the `linx.svs.tsv` file.
#'
#' @param x Path to `linx.svs.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.svs.tsv", package = "gpgr")
#' (l <- linx_svs_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "localTICountEnd")
#'
#' @export
linx_svs_read <- function(x) {
  nm <- c(
    vcfId = "c", svId = "c", clusterId = "d", clusterReason = "c",
    fragileSiteStart = "c", fragileSiteEnd = "c", isFoldback = "c",
    lineTypeStart = "c", lineTypeEnd = "c", junctionCopyNumberMin = "d",
    junctionCopyNumberMax = "d", geneStart = "c", geneEnd = "c",
    replicationTimingStart = "d", replicationTimingEnd = "d",
    localTopologyIdStart = "d", localTopologyIdEnd = "d",
    localTopologyStart = "c", localTopologyEnd = "c",
    localTICountStart = "d", localTICountEnd = "d"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d
}

#' Read LINX Breakend File
#'
#' Reads the `linx.breakend.tsv` file.
#'
#' @param x Path to `linx.breakend.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.breakend.tsv", package = "gpgr")
#' (l <- linx_breakend_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "junctionCopyNumber")
#'
#' @export
linx_breakend_read <- function(x) {
  nm <- c(
    id = "c", svId = "c", isStart = "c",
    gene = "c", transcriptId = "c", canonical = "c",
    geneOrientation = "c", disruptive = "c", reportedDisruption = "c",
    undisruptedCopyNumber = "d", regionType = "c", codingContext = "c",
    biotype = "c", exonicBasePhase = "d", nextSpliceExonRank = "d",
    nextSpliceExonPhase = "d", nextSpliceDistance = "d", totalExonCount = "d",
    type = "c", chromosome = "c", orientation = "c",
    strand = "c", chrBand = "c", exonUp = "d",
    exonDown = "d", junctionCopyNumber = "d"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d
}

#' Read LINX Clusters File
#'
#' Reads the `linx.clusters.tsv` file.
#'
#' @param x Path to `linx.clusters.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.clusters.tsv", package = "gpgr")
#' (l <- linx_clusters_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "clusterDesc")
#'
#' @export
linx_clusters_read <- function(x) {
  nm <- c(
    clusterId = "c", category = "c", synthetic = "c", resolvedType = "c",
    clusterCount = "d", clusterDesc = "c"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d
}

#' Read LINX Links File
#'
#' Reads the `linx.links.tsv` file.
#'
#' @param x Path to `linx.links.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.links.tsv", package = "gpgr")
#' (l <- linx_links_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "ecDna")
#'
#' @export
linx_links_read <- function(x) {
  nm <- c(
    clusterId = "c", chainId = "c", chainIndex = "c", chainCount = "d",
    lowerBreakendId = "d", upperBreakendId = "d", lowerBreakendIsStart = "c",
    upperBreakendIsStart = "c", chromosome = "c", arm = "c", assembled = "c",
    traversedSVCount = "d", length = "d", junctionCopyNumber = "d",
    junctionCopyNumberUncertainty = "d", pseudogeneInfo = "c", ecDna = "c"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d
}

#' Read LINX Fusion File
#'
#' Reads the `linx.fusion.tsv` file.
#'
#' @param x Path to `linx.fusion.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.fusion.tsv", package = "gpgr")
#' (l <- linx_fusion_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "junctionCopyNumber")
#'
#' @export
linx_fusion_read <- function(x) {
  nm <- c(
    fivePrimeBreakendId = "d", threePrimeBreakendId = "d", name = "c",
    reported = "c", reportedType = "c", phased = "c", likelihood = "c",
    chainLength = "d", chainLinks = "d", chainTerminated = "c",
    domainsKept = "c", domainsLost = "c", skippedExonsUp = "d",
    skippedExonsDown = "d", fusedExonUp = "d", fusedExonDown = "d",
    geneStart = "c", geneContextStart = "c", transcriptStart = "c",
    geneEnd = "c", geneContextEnd = "c", transcriptEnd = "c",
    junctionCopyNumber = "d"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d
}

#' Read LINX Driver Catalog File
#'
#' Reads the `linx.driver.catalog.tsv` file.
#'
#' @param x Path to `linx.driver.catalog.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.driver.catalog.tsv", package = "gpgr")
#' (l <- linx_drivercatalog_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "maxCopyNumber")
#'
#' @export
linx_drivercatalog_read <- function(x) {
  nm <- c(
    chromosome = "c", chromosomeBand = "c", gene = "c", driver = "c",
    category = "c", likelihoodMethod = "c", driverLikelihood = "d", `NA` = "c",
    missense = "d", nonsense = "d", splice = "d", inframe = "d",
    frameshift = "d", biallelic = "c", minCopyNumber = "d", maxCopyNumber = "d"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d
}

#' Read LINX Drivers File
#'
#' Reads the `linx.drivers.tsv` file.
#'
#' @param x Path to `linx.drivers.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.drivers.tsv", package = "gpgr")
#' (l <- linx_drivers_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "eventType")
#'
#' @export
linx_drivers_read <- function(x) {
  nm <- c(clusterId = "c", gene = "c", eventType = "c")
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d
}


#### 'vis' outputs ####

#' Read LINX VisCopyNumber File
#'
#' Reads the `linx.vis_copy_number.tsv` file.
#'
#' @param x Path to `linx.vis_copy_number.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_copy_number.tsv", package = "gpgr")
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
  d
}

#' Read LINX VisFusion File
#'
#' Reads the `linx.vis_fusion.tsv` file.
#'
#' @param x Path to `linx.vis_fusion.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_fusion.tsv", package = "gpgr")
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
#' x <- system.file("extdata/linx/tables/linx.vis_gene_exon.tsv", package = "gpgr")
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
  d
}

#' Read LINX VisProteinDomain File
#'
#' Reads the `linx.vis_protein_domain.tsv` file.
#'
#' @param x Path to `linx.vis_protein_domain.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_protein_domain.tsv", package = "gpgr")
#' (l <- linx_visproteindomain_read(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "Info")
#'
#' @export
linx_visproteindomain_read <- function(x) {
  nm <- c(
    SampleId = "c", ClusterId = "c", Transcript = "c", Chromosome = "c",
    Start = "i", End = "i", Info = "c"
  )
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d
}

#' Read LINX VisSegments File
#'
#' Reads the `linx.vis_segments.tsv` file.
#'
#' @param x Path to `linx.vis_segments.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_segments.tsv", package = "gpgr")
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
  d
}

#' Read LINX VisSvData File
#'
#' Reads the `linx.vis_sv_data.tsv` file.
#'
#' @param x Path to `linx.vis_sv_data.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.vis_sv_data.tsv", package = "gpgr")
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
  d
}
