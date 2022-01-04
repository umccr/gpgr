LINX_DESCRIPTIONS <- tibble::tribble(
  ~Table, ~Field, ~Type, ~Description,
  "svs", "vcfId", "c", "ID of break junction for mapping to GRIDSS / PURPLE vcf",
  "svs", "svId", "c", "ID of break junction",
  "svs", "clusterId", "c", "ID of cluster which break junction is assigned to",
  "svs", "clusterReason", "c", paste0(
    "Reason for clustering and svId of clustered break ",
    "junction for other break junction(s) to which the variant has been clustered"
  ),
  "svs", "fragileSiteStart", "c", "Start breakend of break junction is in a known fragile site (T/F)",
  "svs", "fragileSiteEnd", "c", "End breakend of break junction is in a known fragile site (T/F)",
  "svs", "isFoldback", "c", "Break junction is classified as a foldback (T/F)",
  "svs", "lineTypeStart", "c", "Start breakend of break junction is in a known or suspected line source region",
  "svs", "lineTypeEnd", "c", "End breakend of break junction is in a known or suspected line source region",
  "svs", "junctionCopyNumberMin", "d", "Minimum bound JCN estimate for breakjunction",
  "svs", "junctionCopyNumberMax", "d", "Maximum bound JCN estimate for breakjunction",
  "svs", "geneStart", "c", "Gene(s) overlapping start breakend of SV",
  "svs", "geneEnd", "c", "Gene(s) overlapping end breakend of SV",
  "svs", "replicationTimingStart", "d", "-",
  "svs", "replicationTimingEnd", "d", "-",
  "svs", "localTopologyIdStart", "c", paste0(
    "ID for group of proximate breakends to the ",
    "start breakend of break junction within an extending 5kb window"
  ),
  "svs", "localTopologyIdEnd", "c", paste0(
    "ID for group of proximate breakends to the end ",
    "breakend of break junction within an extending 5kb window"
  ),
  "svs", "localTopologyStart", "c", paste0(
    "Local breakend topology type at site of start breakend. ",
    "[ISOLATED_BE, DSB, TI_ONLY, SIMPLE_DUP, FOLDBACK, FOLDBACK_DSB, ",
    "SAME_ORIENT, COMPLEX_FOLDBACK, COMPLEX_LINE, COMPLEX_OTHER]"
  ),
  "svs", "localTopologyEnd", "c", paste0(
    "Local breakend topology type at site of end breakend. ",
    "[ISOLATED_BE, DSB, TI_ONLY, SIMPLE_DUP, FOLDBACK, FOLDBACK_DSB, ",
    "SAME_ORIENT, COMPLEX_FOLDBACK, COMPLEX_LINE, COMPLEX_OTHER]"
  ),
  "svs", "localTICountStart", "d", "Number of chained templated insertions in local topology group of start breakend",
  "svs", "localTICountEnd", "d", "Number of chained templated insertions in local topology group of end breakend",
  "breakend", "id", "c", "Id of breakend annotation",
  "breakend", "svId", "c", "Id of break junction",
  "breakend", "isStart", "c", "Annotation relates to the start breakend of the break junction (1 = true,0 = false)",
  "breakend", "gene", "c", "Gene annotated",
  "breakend", "transcriptId", "c", "Ensembl stable transcript id of annotation",
  "breakend", "canonical", "c", paste0(
    "Transcript is the canonical transcript of the gene. ",
    "Linx annotates 1 record for each canonical transcript ",
    "overlapping the breakend + a record for any non-canonical ",
    "transcript that is annotated as part of a fusion"
  ),
  "breakend", "geneOrientation", "c", paste0(
    "Orientation which breakend points relative to the ",
    "gene taking into account both gene strand and breakend orientation."
  ),
  "breakend", "disruptive", "c", "Breakend is part of a break junction which disrupts the exonic sequence of the transcript",
  "breakend", "reportedDisruption", "c", "Breakend is disruptive and gene is flagged as reportable for disruption",
  "breakend", "undisruptedCopyNumber", "d", paste0(
    "Number of remaining wildtype alleles of the gene that ",
    "are not disrupted by the breakend.  If <0.5 then disruption is considered Homozygous"
  ),
  "breakend", "regionType", "c", paste0(
    "Location of the breakend relative to the transcript. ",
    "[UPSTREAM (within 10kb upstream of the 1st base of the transcript), INTRONIC, EXONIC]"
  ),
  "breakend", "codingContext", "c", paste0(
    "Location of the breakend relative to the coding context of the transcript. ",
    "[CODING, NON_CODING, UTR_5P, UTR_3P, ENHANCER (IG enhancer rearrangements only)"
  ),
  "breakend", "biotype", "c", "Ensembl biotype of the transcript",
  "breakend", "exonicBasePhase", "d", "If regionType = EXONIC, the exact base phase of the breakend, else -1",
  "breakend", "nextSpliceExonRank", "d", paste0(
    "The exon rank of the next splice acceptor (if gene ",
    "orientation is 'DOWNSTREAM') or splice donor (if gene orientation is 'UPSTREAM')"
  ),
  "breakend", "nextSpliceExonPhase", "d", paste0(
    "The phase of the 1st base after the next splice acceptor ",
    "(if gene orientation is 'DOWNSTREAM') or splice donor (if gene orientation is 'UPSTREAM')"
  ),
  "breakend", "nextSpliceDistance", "d", "The distance in bases to the next splice site identified in nextSpliceExonRank",
  "breakend", "totalExonCount", "d", "Total number of exons in the transcript",
  "breakend", "type", "c", "-",
  "breakend", "chromosome", "c", "-",
  "breakend", "orientation", "c", "-",
  "breakend", "strand", "c", "-",
  "breakend", "chrBand", "c", "-",
  "breakend", "exonUp", "d", "-",
  "breakend", "exonDown", "d", "-",
  "breakend", "junctionCopyNumber", "d", "-",
  "clusters", "clusterId", "c", "Unique Id for the cluster",
  "clusters", "category", "c", "High level categorisation of the cluster classification",
  "clusters", "synthetic", "c", paste0(
    "Set to TRUE if the cluster is resolved to a non complex ",
    "type by simplification of a short templated insertion (<1kb)"
  ),
  "clusters", "resolvedType", "c", "Resolved classification of the cluster.",
  "clusters", "clusterCount", "d", "The number of break junctions in the cluster",
  "clusters", "clusterDesc", "c", "String containing the types and counts of break junctions in the cluster. eg. DEL=2_INV=2",
  "links", "clusterId", "c", "Id of the cluster which contains the link",
  "links", "chainId", "c", paste0(
    "Id of the chain to which the link belongs representing a ",
    "multi-segment prediction of the derivative chromosome"
  ),
  "links", "chainIndex", "c", paste0(
    "Position of the linked segment in the chain. The predicted ",
    "chain can be reconstructed by traversing each linked segment in order ie. 0,1,...,n"
  ),
  "links", "chainCount", "d", "Total count of linked segments in the chan",
  "links", "lowerBreakendId", "c", "svId of the leftmost breakend of the linked segment",
  "links", "upperBreakendId", "c", "svId of the rightmost breakend of the linked segment",
  "links", "lowerBreakendIsStart", "c", "True if the lower breakend is the start breakend of the break junction",
  "links", "upperBreakendIsStart", "c", "True if the right breakend is the start breakend of the break junction",
  "links", "chromosome", "c", "Chromosome of the linked segment",
  "links", "arm", "c", "Arm (P/Q) of the linked segment",
  "links", "assembled", "c", "True if the segment is linked by a GRIDSS assembly",
  "links", "traversedSVCount", "d", "The number of other breakends that are located on the linked segment",
  "links", "length", "d", "Length of the linked segment",
  "links", "junctionCopyNumber", "d", "Predicted copy number of the chain",
  "links", "junctionCopyNumberUncertainty", "d", "Uncertainty in the copy number of the chain",
  "links", "pseudogeneInfo", "c", paste0(
    "If the segment precisely matches an exon of an ensembl gene, ",
    "then contains details of the matching exon:  {geneName;TranscriptId,ExonRank,ExonLength}"
  ),
  "links", "ecDna", "c", "True if the link is predicted to be part of a DM / ecDNA chain",
  "fusion", "fivePrimeBreakendId", "c", "Id of the 5' breakend in the fusion",
  "fusion", "threePrimeBreakendId", "c", "Id of the 3' breakend in the fusion",
  "fusion", "name", "c", "Name of the fusion in the form 5'GENE_3'GENE",
  "fusion", "reported", "c", "True if the fusion meets all reportable fusion criteria for Linx",
  "fusion", "reportedType", "c", paste0(
    "If one or both of the genes matches  a promiscuous gene or known rearrangement ",
    "in the HMF fusion knowledgebase, then the type of reportable gene pair: ",
    "[KNOWN_PAIR, PROMISCUOUS_5, PROMISCUOUS_3, PROMISCUOUS_BOTH, EXON_DEL_DUP, ",
    "IG_PROMISCUOUS, IG_KNOWN_PAIR, KNOWN_PAIR_UNMMABLE_3 or NONE (if no match is found)"
  ),
  "fusion", "phased", "c", "Set to 1 if a phased fusion can be formed (after allowing for exon skipping)",
  "fusion", "likelihood", "c", "-",
  "fusion", "chainLength", "d", "0 for simple fusions.  If fusion is chained equal to the total length of segments chained between 5' and 3' partners",
  "fusion", "chainLinks", "d", "0 for simple fusions.  If fusion is chained equal to the number of segments chained between 5' and 3' partners",
  "fusion", "chainTerminated", "c", paste0(
    "True if the fusion is interrupted either on the 5’ partner ",
    "side by a chained breakend prior to the start of the 5’ gene ",
    "or by a chained breakend prior to the last coding base of the 3’ gene"
  ),
  "fusion", "domainsKept", "c", "List of 3' partner domains retained in fusion product (as annotated by PROSITE profiles)",
  "fusion", "domainsLost", "c", "List of 3' partner domains lost in fusion product (as annotated by PROSITE profiles)",
  "fusion", "skippedExonsUp", "d", "Count of splice donors required to be skipped on 5' partner side to form an inframe fusion.",
  "fusion", "skippedExonsDown", "d", "Count of splice donors required to be skipped on 3' partner side to form an inframe fusion",
  "fusion", "fusedExonUp", "d", "Last exon fused on 5' partner side",
  "fusion", "fusedExonDown", "d", "First exon fused on 3' partner side",
  "fusion", "geneStart", "c", "-",
  "fusion", "geneContextStart", "c", "-",
  "fusion", "transcriptStart", "c", "-",
  "fusion", "geneEnd", "c", "-",
  "fusion", "geneContextEnd", "c", "-",
  "fusion", "transcriptEnd", "c", "-",
  "fusion", "junctionCopyNumber", "d", "-",
  "drivercatalog", "chromosome", "c", "Chromosome of gene",
  "drivercatalog", "chromosomeBand", "c", "Chromosome band of driver",
  "drivercatalog", "gene", "c", "Gene name",
  "drivercatalog", "driver", "c", "Driver type [AMP, DEL, MUTATION]",
  "drivercatalog", "category", "c", "Gene driver type [ONCO, TSG]",
  "drivercatalog", "likelihoodMethod", "c", "Method used to determine likelihood [AMP, DEL, BIALLELIC, DNDS, HOTSPOT, INFRAME]",
  "drivercatalog", "driverLikelihood", "d", "Likelihood that gene is a driver",
  "drivercatalog", "NA", "c", "-",
  "drivercatalog", "missense", "d", "Number of missense variants in gene",
  "drivercatalog", "nonsense", "d", "Number of nonsense variants in gene",
  "drivercatalog", "splice", "d", "Number of splice variants in gene",
  "drivercatalog", "inframe", "d", "Number of inframe variants in gene",
  "drivercatalog", "frameshift", "d", "Number of frameshift variants in gene",
  "drivercatalog", "biallelic", "c", "True if any variants in the gene are biallelic",
  "drivercatalog", "minCopyNumber", "d", "Minimum copy number found in the gene exons",
  "drivercatalog", "maxCopyNumber", "d", "Maximum copy number found in the gene exons",
  "drivers", "clusterId", "c", "Id of cluster which break junction associated with driver. Set to -1 for ARM or CHR level events.",
  "drivers", "gene", "c", "Gene of driver. Multiple clusters may be linked to a gene for a sample",
  "drivers", "eventType", "c", paste0(
    "Type of driver. [GAIN (amplification by SV), GAIN_ARM (amplification of whole arm), ",
    "GAIN_CHR (amplification of whole chromosome), DEL (homozygous deletion), ",
    "LOH (focal LOH), LOH_ARM (arm level LOH), LOH_CHR (chromosome level LOH), ",
    "LOH_SV_TELO (LOH from SV to telomere), LOH_SV_CENTRO (LOH from SV to centromere), ",
    "HOM_DUP_DISRUPTION (homozygous disruption via cross exonic tandem duplication), ",
    "HOM_DEL_DISRUPTION (homozygous disruption without homozygous copy number loss)]"
  )
)

#' Read LINX Table
#'
#' Generic function to read LINX tables.
#'
#' @param x Path to LINX file.
#' @param tab Table type to subset on.
#' @keywords internal
#' @noRd
linx_read_table <- function(x, tab) {
  subset_tab <- function(tab) {
    LINX_DESCRIPTIONS |>
      dplyr::filter(.data$Table == tab) |>
      dplyr::select(-.data$Table)
  }

  descr <- subset_tab(tab)
  nm <- descr[["Type"]]
  cnames <- descr[["Field"]]
  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == cnames))
  d
}

#' Get LINX Description for Table Type
#'
#' Gets LINX Description for Table Type
#'
#' @param tab Table type to subset on.
#' @keywords internal
#' @noRd
linx_descr_tab <- function(tab) {
  LINX_DESCRIPTIONS |>
    dplyr::filter(.data$Table == tab) |>
    dplyr::select(.data$Field, .data$Description)
}

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
  linx_read_table(x, "svs")
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
  linx_read_table(x, "breakend")
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
  linx_read_table(x, "clusters")
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
  linx_read_table(x, "links")
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
  linx_read_table(x, "fusion")
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
  linx_read_table(x, "drivercatalog")
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
  linx_read_table(x, "drivers")
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
  descr <- linx_descr_tab("svs")
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
  descr <- linx_descr_tab("breakend")
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
  descr <- linx_descr_tab("clusters")
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
  descr <- linx_descr_tab("links")
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
  descr <- linx_descr_tab("fusion")
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
  descr <- linx_descr_tab("drivercatalog")
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
  descr <- linx_descr_tab("drivers")
  list(tab = d, descr = descr)
}
