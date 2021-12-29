#' Read LINX SVS File
#'
#' Reads the `linx.svs.tsv` file.
#'
#' @param x Path to `linx.svs.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.svs.tsv.gz", package = "gpgr")
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
#' x <- system.file("extdata/linx/tables/linx.breakend.tsv.gz", package = "gpgr")
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
#' x <- system.file("extdata/linx/tables/linx.clusters.tsv.gz", package = "gpgr")
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
#' x <- system.file("extdata/linx/tables/linx.links.tsv.gz", package = "gpgr")
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
#' x <- system.file("extdata/linx/tables/linx.fusion.tsv.gz", package = "gpgr")
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
#' x <- system.file("extdata/linx/tables/linx.driver.catalog.tsv.gz", package = "gpgr")
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
#' x <- system.file("extdata/linx/tables/linx.drivers.tsv.gz", package = "gpgr")
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

#---- 'vis' outputs ----#

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
#' x <- system.file("extdata/linx/tables/linx.vis_protein_domain.tsv.gz", package = "gpgr")
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
  d
}

#' Process LINX SVS File
#'
#' Processes the `linx.svs.tsv` file.
#'
#' @param x Path to `linx.svs.tsv` file.
#'
#' @return The input file as a processed tibble.
#' @export
linx_svs_process <- function(x) {
  d <- linx_svs_read(x)
  descr <- tibble::tribble(
    ~Field, ~Description,
    "vcfId", "ID of break junction for mapping to GRIDSS / PURPLE vcf",
    "svId", "ID of break junction",
    "clusterId", "ID of cluster which break junction is assigned to",
    "clusterReason", paste0(
      "Reason for clustering and svId of clustered break ",
      "junction for other break junction(s) to which the variant has been clustered"
    ),
    "fragileSiteStart", "Start breakend of break junction is in a known fragile site (T/F)",
    "fragileSiteEnd", "End breakend of break junction is in a known fragile site (T/F)",
    "isFoldback", "Break junction is classified as a foldback (T/F)",
    "lineTypeStart", "Start breakend of break junction is in a known or suspected line source region",
    "lineTypeEnd", "End breakend of break junction is in a known or suspected line source region",
    "junctionCopyNumberMin", "Minimum bound JCN estimate for breakjunction",
    "junctionCopyNumberMax", "Maximum bound JCN estimate for breakjunction",
    "geneStart", "Gene(s) overlapping start breakend of SV",
    "geneEnd", "Gene(s) overlapping end breakend of SV",
    "replicationTimingStart", "-",
    "replicationTimingEnd", "-",
    "localTopologyIdStart", paste0(
      "ID for group of proximate breakends to the ",
      "start breakend of break junction within an extending 5kb window"
    ),
    "localTopologyIdEnd", paste0(
      "ID for group of proximate breakends to the end ",
      "breakend of break junction within an extending 5kb window"
    ),
    "localTopologyStart", paste0(
      "Local breakend topology type at site of start breakend. ",
      "[ISOLATED_BE, DSB, TI_ONLY, SIMPLE_DUP, FOLDBACK, FOLDBACK_DSB, ",
      "SAME_ORIENT, COMPLEX_FOLDBACK, COMPLEX_LINE, COMPLEX_OTHER]"
    ),
    "localTopologyEnd", paste0(
      "Local breakend topology type at site of end breakend. ",
      "[ISOLATED_BE, DSB, TI_ONLY, SIMPLE_DUP, FOLDBACK, FOLDBACK_DSB, ",
      "SAME_ORIENT, COMPLEX_FOLDBACK, COMPLEX_LINE, COMPLEX_OTHER]"
    ),
    "localTICountStart", "Number of chained templated insertions in local topology group of start breakend",
    "localTICountEnd", "Number of chained templated insertions in local topology group of end breakend"
  )
  assertthat::assert_that(
    all(colnames(d) %in% descr[["Field"]]),
    ncol(d) == nrow(descr)
  )

  list(tab = d, descr = descr)
}

#' Process LINX Breakend File
#'
#' Processes the `linx.breakend.tsv` file.
#'
#' @param x Path to `linx.breakend.tsv` file.
#'
#' @return The input file as a tibble with column description.
#' @export
linx_breakend_process <- function(x) {
  d <- linx_breakend_read(x)
  descr <- tibble::tribble(
    ~Field, ~Description,
    "id", "Id of breakend annotation",
    "svId", "Id of break junction",
    "isStart", "Annotation relates to the start breakend of the break junction (1 = true,0 = false)",
    "gene", "Gene annotated",
    "transcriptId", "Ensembl stable transcript id of annotation",
    "canonical", paste0(
      "Transcript is the canonical transcript of the gene. ",
      "Linx annotates 1 record for each canonical transcript ",
      "overlapping the breakend + a record for any non-canonical ",
      "transcript that is annotated as part of a fusion"
    ),
    "geneOrientation", paste0(
      "Orientation which breakend points relative to the ",
      "gene taking into account both gene strand and breakend orientation."
    ),
    "disruptive", "Breakend is part of a break junction which disrupts the exonic sequence of the transcript",
    "reportedDisruption", "Breakend is disruptive and gene is flagged as reportable for disruption",
    "undisruptedCopyNumber", paste0(
      "Number of remaining wildtype alleles of the gene that ",
      "are not disrupted by the breakend.  If <0.5 then disruption is considered Homozygous"
    ),
    "regionType", paste0(
      "Location of the breakend relative to the transcript. ",
      "[UPSTREAM (within 10kb upstream of the 1st base of the transcript), INTRONIC, EXONIC]"
    ),
    "codingContext", paste0(
      "Location of the breakend relative to the coding context of the transcript. ",
      "[CODING, NON_CODING, UTR_5P, UTR_3P, ENHANCER (IG enhancer rearrangements only)"
    ),
    "biotype", "Ensembl biotype of the transcript",
    "exonicBasePhase", "If regionType = EXONIC, the exact base phase of the breakend, else -1",
    "nextSpliceExonRank", paste0(
      "The exon rank of the next splice acceptor (if gene ",
      "orientation is 'DOWNSTREAM') or splice donor (if gene orientation is 'UPSTREAM')"
    ),
    "nextSpliceExonPhase", paste0(
      "The phase of the 1st base after the next splice acceptor ",
      "(if gene orientation is 'DOWNSTREAM') or splice donor (if gene orientation is 'UPSTREAM')"
    ),
    "nextSpliceDistance", "The distance in bases to the next splice site identified in nextSpliceExonRank",
    "totalExonCount", "Total number of exons in the transcript",
    "type", "-",
    "chromosome", "-",
    "orientation", "-",
    "strand", "-",
    "chrBand", "-",
    "exonUp", "-",
    "exonDown", "-",
    "junctionCopyNumber", "-"
  )
  assertthat::assert_that(
    all(colnames(d) %in% descr[["Field"]]),
    ncol(d) == nrow(descr)
  )
  list(tab = d, descr = descr)
}

#' Process LINX Clusters File
#'
#' Processes the `linx.clusters.tsv` file.
#'
#' @param x Path to `linx.clusters.tsv` file.
#'
#' @return The input file as a tibble.
#' @export
linx_clusters_process <- function(x) {
  d <- linx_clusters_read(x)
  descr <- tibble::tribble(
    ~Field, ~Description,
    "clusterId", "Unique Id for the cluster",
    "category", "High level categorisation of the cluster classification",
    "synthetic", paste0(
      "Set to TRUE if the cluster is resolved to a non complex ",
      "type by simplification of a short templated insertion (<1kb)"
    ),
    "resolvedType", "Resolved classification of the cluster.",
    "clusterCount", "The number of break junctions in the cluster",
    "clusterDesc", "String containing the types and counts of break junctions in the cluster. eg. DEL=2_INV=2"
  )
  assertthat::assert_that(
    all(colnames(d) %in% descr[["Field"]]),
    ncol(d) == nrow(descr)
  )
  list(tab = d, descr = descr)
}

#' Process LINX Links File
#'
#' Processes the `linx.links.tsv` file.
#'
#' @param x Path to `linx.links.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.links.tsv.gz", package = "gpgr")
#' (l <- linx_links_process(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "ecDna")
#'
#' @export
linx_links_process <- function(x) {
  d <- linx_links_read(x)
  descr <- tibble::tribble(
    ~Field, ~Description,
    "clusterId", "Id of the cluster which contains the link",
    "chainId", paste0(
      "Id of the chain to which the link belongs representing a ",
      "multi-segment prediction of the derivative chromosome"
    ),
    "chainIndex", paste0(
      "Position of the linked segment in the chain. The predicted ",
      "chain can be reconstructed by traversing each linked segment in order ie. 0,1,...,n"
    ),
    "chainCount", "Total count of linked segments in the chan",
    "lowerBreakendId", "svId of the leftmost breakend of the linked segment",
    "upperBreakendId", "svId of the rightmost breakend of the linked segment",
    "lowerBreakendIsStart", "True if the lower breakend is the start breakend of the break junction",
    "upperBreakendIsStart", "True if the right breakend is the start breakend of the break junction",
    "chromosome", "Chromosome of the linked segment",
    "arm", "Arm (P/Q) of the linked segment",
    "assembled", "True if the segment is linked by a GRIDSS assembly",
    "traversedSVCount", "The number of other breakends that are located on the linked segment",
    "length", "Length of the linked segment",
    "junctionCopyNumber", "Predicted copy number of the chain",
    "junctionCopyNumberUncertainty", "Uncertainty in the copy number of the chain",
    "pseudogeneInfo", paste0(
      "If the segment precisely matches an exon of an ensembl gene, ",
      "then contains details of the matching exon:  {geneName;TranscriptId,ExonRank,ExonLength}"
    ),
    "ecDna", "True if the link is predicted to be part of a DM / ecDNA chain"
  )
  assertthat::assert_that(
    all(colnames(d) %in% descr[["Field"]]),
    ncol(d) == nrow(descr)
  )
  list(tab = d, descr = descr)
}

#' Process LINX Fusion File
#'
#' Processes the `linx.fusion.tsv` file.
#'
#' @param x Path to `linx.fusion.tsv` file.
#'
#' @return The input file as a tibble.
#' @export
linx_fusion_process <- function(x) {
  d <- linx_fusion_read(x)
  descr <- tibble::tribble(
    ~Field, ~Description,
    "fivePrimeBreakendId", "Id of the 5' breakend in the fusion",
    "threePrimeBreakendId", "Id of the 3' breakend in the fusion",
    "name", "Name of the fusion in the form 5'GENE_3'GENE",
    "reported", "True if the fusion meets all reportable fusion criteria for Linx",
    "reportedType", paste0(
      "If one or both of the genes matches  a promiscuous gene or known rearrangement ",
      "in the HMF fusion knowledgebase, then the type of reportable gene pair: ",
      "[KNOWN_PAIR, PROMISCUOUS_5, PROMISCUOUS_3, PROMISCUOUS_BOTH, EXON_DEL_DUP, ",
      "IG_PROMISCUOUS, IG_KNOWN_PAIR, KNOWN_PAIR_UNMMABLE_3 or NONE (if no match is found)"
    ),
    "phased", "Set to 1 if a phased fusion can be formed (after allowing for exon skipping)",
    "likelihood", "-",
    "chainLength", "0 for simple fusions.  If fusion is chained equal to the total length of segments chained between 5' and 3' partners",
    "chainLinks", "0 for simple fusions.  If fusion is chained equal to the number of segments chained between 5' and 3' partners",
    "chainTerminated", paste0(
      "True if the fusion is interrupted either on the 5’ partner ",
      "side by a chained breakend prior to the start of the 5’ gene ",
      "or by a chained breakend prior to the last coding base of the 3’ gene"
    ),
    "domainsKept", "List of 3' partner domains retained in fusion product (as annotated by PROSITE profiles)",
    "domainsLost", "List of 3' partner domains lost in fusion product (as annotated by PROSITE profiles)",
    "skippedExonsUp", "Count of splice donors required to be skipped on 5' partner side to form an inframe fusion.",
    "skippedExonsDown", "Count of splice donors required to be skipped on 3' partner side to form an inframe fusion",
    "fusedExonUp", "Last exon fused on 5' partner side",
    "fusedExonDown", "First exon fused on 3' partner side",
    "geneStart", "-",
    "geneContextStart", "-",
    "transcriptStart", "-",
    "geneEnd", "-",
    "geneContextEnd", "-",
    "transcriptEnd", "-",
    "junctionCopyNumber", "-"
  )
  assertthat::assert_that(
    all(colnames(d) %in% descr[["Field"]]),
    ncol(d) == nrow(descr)
  )
  list(tab = d, descr = descr)
}

#' Process LINX Driver Catalog File
#'
#' Processes the `linx.driver.catalog.tsv` file.
#'
#' @param x Path to `linx.driver.catalog.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.driver.catalog.tsv.gz", package = "gpgr")
#' (l <- linx_drivercatalog_process(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "maxCopyNumber")
#'
#' @export
linx_drivercatalog_process <- function(x) {
  d <- linx_drivercatalog_read(x)
  descr <- tibble::tribble(
    ~Field, ~Description,
    "chromosome", "Chromosome of gene",
    "chromosomeBand", "Chromosome band of driver",
    "gene", "Gene name",
    "driver", "Driver type [AMP, DEL, MUTATION]",
    "category", "Gene driver type [ONCO, TSG]",
    "likelihoodMethod", "Method used to determine likelihood [AMP, DEL, BIALLELIC, DNDS, HOTSPOT, INFRAME]",
    "driverLikelihood", "Likelihood that gene is a driver",
    "NA", "-",
    "missense", "Number of missense variants in gene",
    "nonsense", "Number of nonsense variants in gene",
    "splice", "Number of splice variants in gene",
    "inframe", "Number of inframe variants in gene",
    "frameshift", "Number of frameshift variants in gene",
    "biallelic", "True if any variants in the gene are biallelic",
    "minCopyNumber", "Minimum copy number found in the gene exons",
    "maxCopyNumber", "Maximum copy number found in the gene exons"
  )
  assertthat::assert_that(
    all(colnames(d) %in% descr[["Field"]]),
    ncol(d) == nrow(descr)
  )
  list(tab = d, descr = descr)
}

#' Process LINX Drivers File
#'
#' Processes the `linx.drivers.tsv` file.
#'
#' @param x Path to `linx.drivers.tsv` file.
#'
#' @return The input file as a tibble.
#' @examples
#' x <- system.file("extdata/linx/tables/linx.drivers.tsv.gz", package = "gpgr")
#' (l <- linx_drivers_process(x))
#' @testexamples
#' expect_equal(colnames(l)[ncol(l)], "eventType")
#'
#' @export
linx_drivers_process <- function(x) {
  d <- linx_drivers_read(x)
  descr <- tibble::tribble(
    ~Field, ~Description,
    "clusterId", "Id of cluster which break junction associated with driver. Set to -1 for ARM or CHR level events.",
    "gene", "Gene of driver. Multiple clusters may be linked to a gene for a sample",
    "eventType", paste0(
      "Type of driver. [GAIN (amplification by SV), GAIN_ARM (amplification of whole arm), ",
      "GAIN_CHR (amplification of whole chromosome), DEL (homozygous deletion), ",
      "LOH (focal LOH), LOH_ARM (arm level LOH), LOH_CHR (chromosome level LOH), ",
      "LOH_SV_TELO (LOH from SV to telomere), LOH_SV_CENTRO (LOH from SV to centromere), ",
      "HOM_DUP_DISRUPTION (homozygous disruption via cross exonic tandem duplication), ",
      "HOM_DEL_DISRUPTION (homozygous disruption without homozygous copy number loss)]"
    )
  )
  assertthat::assert_that(
    all(colnames(d) %in% descr[["Field"]]),
    ncol(d) == nrow(descr)
  )
  list(tab = d, descr = descr)
}
