#' Read PURPLE CNV Gene File
#'
#' Reads the `purple.cnv.gene.tsv` file, which summarises copy number
#' alterations of each gene in the HMF panel
#' (see [this table](https://github.com/hartwigmedical/hmftools/tree/master/purple#gene-copy-number-file)).
#'
#' @param x Path to `purple.cnv.gene.tsv` file.
#'
#' @return The input file as a tibble.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.gene.tsv", package = "gpgr")
#' (p <- purple_cnv_som_gene_read(x))
#' @testexamples
#' expect_equal(colnames(p)[ncol(p)], "minMinorAlleleCopyNumber")
#'
#' @export
purple_cnv_som_gene_read <- function(x) {
  nm <- c(
    "chromosome" = "c", "start" = "i", "end" = "i", "gene" = "c",
    "minCopyNumber" = "d", "maxCopyNumber" = "d", "somaticRegions" = "d",
    "transcriptId" = "c", "isCanonical" = "c", "chromosomeBand" = "c",
    "minRegions" = "d", "minRegionStart" = "i", "minRegionEnd" = "i",
    "minRegionStartSupport" = "c", "minRegionEndSupport" = "c",
    "minRegionMethod" = "c", "minMinorAlleleCopyNumber" = "d",
    "depthWindowCount" = "i"
  )

  ctypes <- paste(nm, collapse = "")
  purple_cnv_gene <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(purple_cnv_gene) == length(nm))
  assertthat::assert_that(all(colnames(purple_cnv_gene) == names(nm)))
  purple_cnv_gene
}

#' Process PURPLE CNV Gene File for UMCCRISE
#'
#' Processes the `purple.cnv.gene.tsv` file. Keeps genes that are in the
#' [UMCCR cancer gene list](https://github.com/umccr/genes/blob/893a655801ce92715f05517b5052e4e81904e870/panels/umccr_2019-03-20.tsv)
#' and selects columns of interest.
#'
#' @param x Path to `purple.cnv.gene.tsv` file.
#' @param g Path to gene file containing at least three columns:
#' * `symbol`: gene name (character).
#' * `tumorsuppressor`: is this gene a tumor suppressor (TRUE/FALSE).
#' * `oncogene`: is this gene an oncogene (TRUE/FALSE).
#'
#' @return List with two elements:
#' * `tab`: Tibble filtered to genes found in  `g`.
#' * `descr`: Description of tibble columns.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.gene.tsv", package = "gpgr")
#' g <- system.file("extdata/ref/umccr_cancer_genes_2019-03-20.tsv", package = "gpgr")
#' (pp <- purple_cnv_som_gene_process(x, g))
#' @testexamples
#' expect_equal(colnames(pp$tab)[ncol(pp$tab)], "minRegSupportStartEndMethod")
#'
#' @export
purple_cnv_som_gene_process <- function(x, g = NULL) {
  purple_cnv_gene <- purple_cnv_som_gene_read(x)
  if (is.null(g)) {
    g <- system.file("extdata/ref/umccr_cancer_genes_2019-03-20.tsv", package = "gpgr")
  }
  genes <-
    readr::read_tsv(g, col_types = readr::cols(
      symbol = "c", oncogene = "l", tumorsuppressor = "l"
    )) |>
    dplyr::select("symbol", "oncogene", "tumorsuppressor")
  oncogenes <- genes |>
    dplyr::filter(.data$oncogene) |>
    dplyr::pull(.data$symbol)
  tsgenes <- genes |>
    dplyr::filter(.data$tumorsuppressor) |>
    dplyr::pull(.data$symbol)

  purple_cnv_gene <- purple_cnv_gene |>
    dplyr::filter(.data$gene %in% genes$symbol) |>
    dplyr::mutate(
      chromosome = as.factor(.data$chromosome),
      transcriptID = paste0(.data$transcriptId),
      minRegStartEnd = paste0(.data$minRegionStart, "-", .data$minRegionEnd),
      minRegSupportStartEndMethod = paste0(
        .data$minRegionStartSupport, "-", .data$minRegionEndSupport,
        " (", .data$minRegionMethod, ")"
      ),
      oncogene = .data$gene %in% oncogenes,
      tsgene = .data$gene %in% tsgenes,
      onco_or_ts = dplyr::case_when(
        .data$oncogene & .data$tsgene ~ "onco+ts",
        .data$oncogene ~ "oncogene",
        .data$tsgene ~ "tsgene",
        TRUE ~ ""
      )
    ) |>
    dplyr::select("gene",
      minCN = "minCopyNumber", maxCN = "maxCopyNumber",
      chrom = "chromosome", "start", "end",
      chrBand = "chromosomeBand", "onco_or_ts",
      "transcriptID", minMinorAlleleCN = "minMinorAlleleCopyNumber",
      somReg = "somaticRegions", minReg = "minRegions",
      "minRegStartEnd", "minRegSupportStartEndMethod"
    )

  descr <- dplyr::tribble(
    ~Column, ~Description,
    "gene", "Name of gene",
    "minCN/maxCN", "Min/Max copy number found in gene exons",
    "chrom/start/end", "Chromosome/start/end location of gene transcript",
    "chrBand", "Chromosome band of the gene",
    "onco_or_ts", "oncogene ('oncogene'), tumor suppressor ('tsgene'), or both ('onco+ts'), as reported by [Cancermine](https://github.com/jakelever/cancermine)",
    "transcriptID", "Ensembl transcript ID (dot version)",
    "minMinorAlleleCN", "Minimum allele ploidy found over the gene exons - useful for identifying LOH events",
    "somReg (somaticRegions)", "Count of somatic copy number regions this gene spans",
    "minReg (minRegions)", "Number of somatic regions inside the gene that share the min copy number",
    "minRegStartEnd", "Start/End base of the copy number region overlapping the gene with the minimum copy number",
    "minRegSupportStartEndMethod", "Start/end support of the CN region overlapping the gene with the min CN (plus determination method)"
  )

  list(
    tab = purple_cnv_gene,
    descr = descr
  )
}

sash_read_cnv_tsv <- function(x) {

  nm <- c(
    "chromosome" = "c",
    "start" = "i",
    "end" = "i",
    "svtype" = "c",
    "baf" = "d",
    "bafCount" = "d",
    "copyNumber" = "d",
    "depthWindowCount" = "i",
    "gcContent" = "d",
    "majorAlleleCopyNumber" = "d",
    "method" = "c",
    "minorAlleleCopyNumber" = "d",
    "segmentEndSupport" = "c",
    "segmentStartSupport" = "c",
    "sv_top_tier" = "i",
    "simple_ann" = "c"
  )

  ctypes <- paste(nm, collapse = "")
  d <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(d) == length(nm))
  assertthat::assert_that(all(colnames(d) == names(nm)))
  d
}

split_chromosome_annotations <- function(x) {
  x.grouped <- x |>
    dplyr::group_by(
      record_type=ifelse(grepl('chrom_[0-9XY]+', Detail), 'chromosome', 'other')
    )

  keys <- x.grouped |>
    dplyr::group_keys() |>
    dplyr::pull(record_type)

  x.split <- x.grouped |>
    dplyr::group_split(.keep=FALSE) |>
    purrr::set_names(keys)

  chromosome <- purrr::pluck(x.split, 'chromosome') |>
    dplyr::select(
      "Event ID" = "event_id",
      "Tier (top)",
      "Start" = "start",
      "End" = "end",
      "SV Type" = "svtype",
      "Effect",
      "Chromosome" = "Detail",
      "PURPLE CN" = "copyNumber",
      "PURPLE CN Min+Maj",
    )

  list(
    chromosome=chromosome,
    other=purrr::pluck(x.split, 'other')
  )
}

split_fusion_annotations <- function(x) {

  fusion_annotations <- c(
    'BidFusG',
    'FusG',
    'bidirectional_gene_fusion',
    'gene_fusion'
  )

  x.grouped <- x |>
    dplyr::group_by(
      record_type=ifelse(Effect %in% fusion_annotations, 'fusion', 'other')
    )

  keys <- x.grouped |>
    dplyr::group_keys() |>
    dplyr::pull(record_type)

  x.split <- x.grouped |>
    dplyr::group_split(.keep=FALSE) |>
    purrr::set_names(keys)

  fusions <- purrr::pluck(x.split, 'fusion') |>
    dplyr::select(
      "Event ID" = "event_id",
      "Annotation ID" = "annotation_id",
      "Tier (top)",
      "Start" = "start",
      "End" = "end",
      "SV Type" = "svtype",
      "Effect",
      "Genes",
      "Transcripts" = "Transcript",
      "Detail",
      "PURPLE CN" = "copyNumber",
      "PURPLE CN Min+Maj",
    )

  list(
    fusions=fusions,
    other=purrr::pluck(x.split, 'other')
  )
}

count_event_genes <- function(.data) {
  .data |>
    dplyr::pull(Genes) |>
    strsplit('&') |>
    unlist() |>
    unique() |>
    length()
}

set_event_type <- function(x) {
  x |>
    dplyr::group_by(event_id) |>
    dplyr::mutate(
      `Gene count (event)` = count_event_genes(dplyr::pick(Genes)),
      event_category = dplyr::case_when(
        `Gene count (event)` == 0 ~ 'none',
        `Gene count (event)` == 1 ~ 'small',
        `Gene count (event)` >= 2 & `Gene count (event)` <= 3 ~ 'medium',
        .default = 'large',
      )
    ) |>
    dplyr::ungroup()
}

set_many_transcripts_cnv <- function(x) {
  # Set many transcripts
  x.tmp <- x |>
    dplyr::rowwise() |>
    dplyr::mutate(
      transcript_count = strsplit(Transcript, '&') |> unlist() |> unique() |> length()
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      many_transcripts = ifelse(transcript_count > 2, 'many_transcripts', 'few_transcripts'),
      many_transcripts = ifelse(event_category %in% c('small', 'medium'), many_transcripts, NA),
    ) |>
    # Sort rows
    dplyr::arrange(`Tier (top)`, Genes, Effect)

  # Build the many transcripts table
  mt <- x.tmp |>
    dplyr::filter(many_transcripts == 'many_transcripts') |>
    # Remove unneeded columns and rename others
    dplyr::select(
      "Event ID" = "event_id",
      "Annotation ID" = "annotation_id",
      "Tier (top)",
      "Start" = "start",
      "End" = "end",
      "Transcript count" = "transcript_count",
      "Transcripts" = "Transcript",
    )

  # Clear transcript where appropriate
  x.ready <- x.tmp |>
    dplyr::mutate(
      Transcript = ifelse(
        many_transcripts == 'few_transcripts' | is.na(many_transcripts),
        Transcript,
        paste0('Many transcripts (', transcript_count, ')')
      )
    ) |>
    dplyr::select(c(-many_transcripts, -transcript_count))

  list(
    many_transcripts=mt,
    cnv=x.ready
  )
}

split_priority <- function(x) {
  x.grouped <- x |>
    dplyr::group_by(
      prioritised = ifelse(Tier < 4, 'prioritised', 'unprioritised')
    )

  keys <- x.grouped |>
    dplyr::group_keys() |>
    dplyr::pull(prioritised)

  x.split <- x.grouped |>
    dplyr::group_split(.keep=FALSE) |>
    purrr::set_names(keys)

  list(
    prioritised = purrr::pluck(x.split, 'prioritised'),
    unprioritised = purrr::pluck(x.split, 'unprioritised', .default=x[NULL, ]) |> dplyr::select(-Transcript)
  )
}

split_event_type <- function(x) {

  x.grouped <- x |>
    dplyr::group_by(event_category)

  keys <- x.grouped |>
    dplyr::group_keys() |>
    dplyr::pull(event_category)

  x.split <- x.grouped |>
    dplyr::group_split(.keep=FALSE) |>
    purrr::set_names(keys)

  list(
    none = purrr::pluck(x.split, 'none'),
    small = purrr::pluck(x.split, 'small'),
    medium = purrr::pluck(x.split, 'medium'),
    large = purrr::pluck(x.split, 'large')
  )
}

collapse_annotations <- function(x) {
  x |>
    dplyr::summarise(
      `Gene count (effect)` = length(unique(unlist(strsplit(Genes, '&')))),
      Genes = paste0(sort(unique(unlist(strsplit(Genes, '&')))), collapse=','),
      Detail = paste0(sort(unique(Detail)), collapse=','),
      `Tier (highest)` = min(Tier),
      `Tier (lowest)` = max(Tier),
    ) |>
    dplyr::ungroup()
}

set_many_genes <- function(x) {
  # Set many genes
  x.tmp <- x |>
    dplyr::rowwise() |>
    dplyr::mutate(
      gene_count = strsplit(Genes, c(',|&')) |> unlist() |> unique() |> length()
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      many_genes = ifelse(gene_count > 2, 'many_genes', 'few_genes'),
    ) |>
    # Sort rows
    dplyr::arrange(`Tier (top)`, Genes, Effect)

  # Build the many genes table
  mt <- x.tmp |>
    dplyr::filter(many_genes == 'many_genes') |>
    # Remove unneeded columns and rename others
    dplyr::select(
      "Event ID" = "event_id",
      "Tier (top)",
      "Start" = "start",
      "End" = "end",
      "SV Type" = "svtype",
      "Effect",
      "Gene count (event)",
      "Gene count (effect)",
      "Genes",
    )

  # Clear genes where appropriate
  x.ready <- x.tmp |>
    dplyr::mutate(
      Genes = ifelse(
        many_genes == 'few_genes' | is.na(many_genes),
        Genes,
        paste0('Many genes (', gene_count, ')')
      )
    ) |>
    dplyr::select(c(-many_genes, -gene_count))

  list(
    many_genes=mt,
    cnv=x.ready
  )
}

#' @export
process_cnv_tsv <- function(x) {
  # Read input
  cnv.input <- sash_read_cnv_tsv(x)

  # Prepare input
  cnv.ready <- cnv.input |>
    dplyr::mutate(
      chrom_simple = stringr::str_remove(chromosome, 'chr'),
      start = paste(chrom_simple, base::format(start, big.mark = ",", trim = TRUE), sep=":"),
      end = paste(chrom_simple, base::format(end, big.mark = ",", trim = TRUE), sep=":"),
    ) |>
    dplyr::select(-c(
      chromosome,
      chrom_simple,
    ))

  # Melt annotations
  cnv.tmp <- cnv.ready |>
    dplyr::mutate(event_id = dplyr::row_number()) |>
    # Split into individual annotations
    dplyr::mutate(annotation = strsplit(simple_ann, ",")) |>
    # Convert annotation fields into columns
    tidyr::unnest(annotation) |>
    tidyr::separate(
      annotation, c("Event", "Effect", "Genes", "Transcript", "Detail", "Tier"),
      sep = "\\|", convert = FALSE
    ) |>
    # Create new columns and modify existing ones
    dplyr::mutate(
      copyNumber = as.numeric(copyNumber) |> round(2) %>% sprintf("%.2f", .),
      minorAlleleCopyNumber = as.numeric(minorAlleleCopyNumber) |> round(2) %>% sprintf("%.2f", .),
      majorAlleleCopyNumber = as.numeric(majorAlleleCopyNumber) |> round(2) %>% sprintf("%.2f", .),
      'PURPLE CN Min+Maj' = paste0(minorAlleleCopyNumber, "+", majorAlleleCopyNumber),
      "Tier (top)" = paste0(Tier, " (", sv_top_tier, ")"),
    ) |>
    # Set unique annotation ID
    dplyr::mutate(annotation_id = dplyr::row_number()) |>
    # Remove unused columns
    dplyr::select(-c(
      baf,
      bafCount,
      depthWindowCount,
      Event,
      gcContent,
      majorAlleleCopyNumber,
      method,
      minorAlleleCopyNumber,
      segmentEndSupport,
      segmentStartSupport,
      simple_ann,
      sv_top_tier,
    )) |>
    # Sort rows
    dplyr::arrange(`Tier (top)`, Genes, Effect)

  # Drop specific annotations
  cnv.tmp <- cnv.tmp |>
    dplyr::filter(Effect != 'intergenic_region')

  # Abbreviate effects
  abbreviate_effectv <- Vectorize(gpgr::abbreviate_effect)
  cnv.tmp$Effect <- abbreviate_effectv(cnv.tmp$Effect)

  # Split chromosome annotations
  cnv.chromosome_split <- split_chromosome_annotations(cnv.tmp)
  cnv.chromosome <- cnv.chromosome_split$chromosome
  cnv.tmp <- cnv.chromosome_split$other

  # Split fusion annotations
  cnv.fusion_split <- split_fusion_annotations(cnv.tmp)
  cnv.fusion <- cnv.fusion_split$fusion
  cnv.tmp <- cnv.fusion_split$other

  # Set event type based on the number of genes impacted
  cnv.tmp <- set_event_type(cnv.tmp)

  # Set and create many transcripts table
  cnv.many_transcripts_data <- set_many_transcripts_cnv(cnv.tmp)
  cnv.many_transcripts <- purrr::pluck(cnv.many_transcripts_data, 'many_transcripts')
  cnv.tmp <- purrr::pluck(cnv.many_transcripts_data, 'cnv')

  # Split into prioritised and unprioritised
  cnv.priority_split <- split_priority(cnv.tmp)
  cnv.unprioritised.tmp <- cnv.priority_split$unprioritised
  cnv.prioritised.tmp <- cnv.priority_split$prioritised

  # Split into event type assigned above
  cnv.event_type_split <- split_event_type(cnv.prioritised.tmp)
  cnv.small.tmp <- cnv.event_type_split$small
  cnv.medium.tmp <- cnv.event_type_split$medium
  cnv.large.tmp <- cnv.event_type_split$large

  # Temporary set dplyr summarise warning off
  summarise_info_opt <- getOption('dplyr.summarise.inform')
  options(dplyr.summarise.inform = FALSE)

  # Collapse large event type annotations into single records, drop transcripts
  cnv.large.tmp <- cnv.large.tmp |>
    dplyr::group_by(dplyr::across(c(-Genes, -Transcript, -Tier, -Detail, -annotation_id))) |>
    collapse_annotations()

  cnv.large.many_genes_data <- set_many_genes(cnv.large.tmp)
  cnv.large.many_genes <- cnv.large.many_genes_data$many_genes
  cnv.large.tmp <- cnv.large.many_genes_data$cnv

  # Select then collapse only large event type annotations into single records
  cnv.unprioritised.tmp <- cnv.unprioritised.tmp |>
    dplyr::filter(event_category == 'large') |>
    dplyr::group_by(dplyr::across(c(-Genes, -Tier, -Detail, -annotation_id))) |>
    collapse_annotations() |>
    dplyr::ungroup() |>
    dplyr::bind_rows(
      cnv.unprioritised.tmp |> dplyr::filter(event_category != 'large')
    )
  cnv.unprioritised.many_genes_data <- set_many_genes(cnv.unprioritised.tmp)
  cnv.unprioritised.many_genes <- cnv.unprioritised.many_genes_data$many_genes
  cnv.unprioritised.tmp <- cnv.unprioritised.many_genes_data$cnv

  # Restore dplyr summarise warning setting
  options(dplyr.summarise.inform = summarise_info_opt)

  # Select columns for remaining tables
  columns <- c(
      "Event ID" = "event_id",
      "Tier (top)",
      "Tier lowest (effect)",
      "Tier lowest (effect)",
      "Start" = "start",
      "End" = "end",
      "SV Type" = "svtype",
      "Effect",
      "Gene count (event)",
      "Gene count (effect)",
      "Genes",
      "Transcripts" = "Transcript",
      "Detail",
      "PURPLE CN" = "copyNumber",
      "PURPLE CN Min+Maj"
  )

  cnv.small <- cnv.small.tmp |> dplyr::select(tidyselect::any_of(columns))
  cnv.medium <- cnv.medium.tmp |> dplyr::select(tidyselect::any_of(columns))
  cnv.large <- cnv.large.tmp |> dplyr::select(tidyselect::any_of(columns))
  cnv.unprioritised <- cnv.unprioritised.tmp |> dplyr::select(tidyselect::any_of(columns))

  list(
    small = cnv.small,
    medium = cnv.medium,
    large = cnv.large,
    many_genes_prioritised = cnv.large.many_genes,
    many_transcripts = cnv.many_transcripts,
    unprioritised = cnv.unprioritised,
    many_genes_unprioritised = cnv.unprioritised.many_genes,
    fusion = cnv.fusion,
    chromsome = cnv.chromosome
  )
}

#' Read PURPLE CNV Somatic File
#'
#' Reads the `purple.cnv.somatic.tsv` file, which contains the copy number
#' profile of all (contiguous) segments of the tumor sample
#' (see [this table](https://github.com/hartwigmedical/hmftools/tree/master/purple#copy-number-file)).
#'
#' @param x Path to `purple.cnv.somatic.tsv` file.
#'
#' @return The input file as a tibble.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' (p <- purple_cnv_som_read(x))
#' @testexamples
#' expect_equal(colnames(p)[ncol(p)], "majorAlleleCopyNumber")
#'
#' @export
purple_cnv_som_read <- function(x) {
  nm <- c(
    "chromosome" = "c", "start" = "i", "end" = "i",
    "copyNumber" = "d", "bafCount" = "d", "observedBAF" = "d",
    "baf" = "d", "segmentStartSupport" = "c", "segmentEndSupport" = "c",
    "method" = "c", "depthWindowCount" = "i", "gcContent" = "d",
    "minStart" = "i", "maxStart" = "i", "minorAlleleCopyNumber" = "d",
    "majorAlleleCopyNumber" = "d"
  )
  ctypes <- paste(nm, collapse = "")
  purple_cnv_somatic <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(purple_cnv_somatic) == length(nm))
  assertthat::assert_that(all(colnames(purple_cnv_somatic) == names(nm)))
  purple_cnv_somatic
}

#' Process PURPLE CNV Somatic File for UMCCRISE
#'
#' Processes the `purple.cnv.somatic.tsv` file and selects columns of interest.
#'
#' @param x Path to `purple.cnv.somatic.tsv` file.
#'
#' @return List with two elements:
#' * `tab`: Tibble with more condensed columns.
#' * `descr`: Description of tibble columns.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' (pp <- purple_cnv_som_process(x))
#' @testexamples
#' expect_equal(colnames(pp$tab)[ncol(pp$tab)], "GC (windowCount)")
#'
#' @export
purple_cnv_som_process <- function(x) {
  purple_cnv_somatic <- purple_cnv_som_read(x)
  purple_cnv_somatic <- purple_cnv_somatic |>
    dplyr::mutate(
      Chr = as.factor(.data$chromosome),
      minorAlleleCopyNumber = round(.data$minorAlleleCopyNumber, 1),
      majorAlleleCopyNumber = round(.data$majorAlleleCopyNumber, 1),
      `CopyNumber Min+Maj` = paste0(.data$minorAlleleCopyNumber, "+", .data$majorAlleleCopyNumber),
      copyNumber = round(.data$copyNumber, 1),
      bafAdj = round(.data$baf, 2),
      gcContent = round(.data$gcContent, 2),
      `Start/End SegSupport` = paste0(.data$segmentStartSupport, "-", .data$segmentEndSupport),
      `BAF (count)` = paste0(.data$bafAdj, " (", .data$bafCount, ")"),
      `GC (windowCount)` = paste0(.data$gcContent, " (", .data$depthWindowCount, ")")
    ) |>
    dplyr::select(
      "Chr",
      Start = "start", End = "end", CN = "copyNumber",
      `CN Min+Maj` = "CopyNumber Min+Maj", "Start/End SegSupport",
      Method = "method", "BAF (count)", "GC (windowCount)"
    )


  descr <- dplyr::tribble(
    ~Column, ~Description,
    "Chr/Start/End", "Coordinates of copy number segment",
    "CN", "Fitted absolute copy number of segment adjusted for purity and ploidy",
    "CN Min+Maj", "CopyNumber of minor + major allele adjusted for purity",
    "Start/End SegSupport", paste0(
      "Type of SV support for the CN breakpoint at ",
      "start/end of region. Allowed values: ",
      "CENTROMERE, TELOMERE, INV, DEL, DUP, BND (translocation), ",
      "SGL (single breakend SV support), NONE (no SV support for CN breakpoint), ",
      "MULT (multiple SV support at exact breakpoint)"
    ),
    "Method", paste0(
      "Method used to determine the CN of the region. Allowed values: ",
      "BAF_WEIGHTED (avg of all depth windows for the region), ",
      "STRUCTURAL_VARIANT (inferred using ploidy of flanking SVs), ",
      "LONG_ARM (inferred from the long arm), GERMLINE_AMPLIFICATION ",
      "(inferred using special logic to handle regions of germline amplification)"
    ),
    "BAF (count)", "Tumor BAF after adjusted for purity and ploidy (Count of AMBER baf points covered by this segment)",
    "GC (windowCount)", "Proportion of segment that is G or C (Count of COBALT windows covered by this segment)"
  )

  list(
    tab = purple_cnv_somatic,
    descr = descr
  )
}

#' Read PURPLE QC file
#'
#' Reads the `purple.qc` file.
#'
#' @param x Path to the `purple.qc` file.
#'
#' @return The input file as a tibble and a summarised tibble with a
#' description of each metric.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.qc", package = "gpgr")
#' (q <- purple_qc_read(x))
#' @testexamples
#' expect_true(q$raw[1, "value", drop = TRUE] == "WARN_DELETED_GENES")
#'
#' @export
purple_qc_read <- function(x) {
  purple_qc <-
    readr::read_tsv(x, col_names = c("key", "value"), col_types = "cc") |>
    dplyr::mutate(value = toupper(.data$value))

  nm <- c(
    "QCStatus", "Method", "CopyNumberSegments",
    "UnsupportedCopyNumberSegments", "Purity", "AmberGender",
    "CobaltGender", "DeletedGenes", "Contamination", "GermlineAberrations",
    "AmberMeanDepth"
  )

  assertthat::assert_that(all(purple_qc$key == nm))
  q <- structure(purple_qc$value, names = purple_qc$key)
  # the n column is used for arranging the final summary table rows in the report
  summary <- dplyr::tribble(
    ~n, ~variable, ~value, ~details,
    1, "QC_Status", glue::glue('{q["QCStatus"]}'),
    paste("See 'Description'."),
    13, "Method", glue::glue('{q["Method"]}'),
    glue::glue("Fit method (NORMAL, HIGHLY_DIPLOID, SOMATIC or NO_TUMOR)."),
    14, "CopyNumberSegments",
    glue::glue('{q["CopyNumberSegments"]} (Unsupported: {q["UnsupportedCopyNumberSegments"]})'),
    "# of CN segments.",
    2, "Purity", glue::glue('{q["Purity"]}'), "",
    17, "Gender", glue::glue('Amber: {q["AmberGender"]}; Cobalt: {q["CobaltGender"]}'), "",
    14, "DeletedGenes", glue::glue('{q["DeletedGenes"]}'), "# of homozygously deleted genes.",
    15, "Contamination", glue::glue('{q["Contamination"]}'),
    "Rate of contamination in tumor sample as determined by AMBER.",
    16, "GermlineAberrations", glue::glue('{q["GermlineAberrations"]}'),
    "Can be one or more of: KLINEFELTER, TRISOMY_X/21/13/18/15, XYY, MOSAIC_X.",
    18, "AmberMeanDepth", glue::glue('{q["AmberMeanDepth"]}'),
    "Mean depth as determined by AMBER.",
  )

  list(
    raw = purple_qc,
    summary = summary
  )
}

#' Read PURPLE Purity file
#'
#' Reads the `purple.purity.tsv` file containing a summary of the purity fit.
#'
#' @param x Path to the `purple.purity.tsv` file.
#'
#' @return The input file as a tibble and a summarised tibble with a
#' description of each metric.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.purity.tsv", package = "gpgr")
#' (p <- purple_purity_read(x))
#' @testexamples
#' expect_equal(p$raw[1, "column", drop = TRUE], "purity")
#' expect_equal(p$raw[nrow(p$raw), "column", drop = TRUE], "svTumorMutationalBurden")
#'
#' @export
purple_purity_read <- function(x) {
  tab <- dplyr::tribble(
    ~column, ~type,
    "purity", "d",
    "normFactor", "d",
    "score", "d",
    "diploidProportion", "d",
    "ploidy", "d",
    "gender", "c",
    "status", "c",
    "polyclonalProportion", "d",
    "minPurity", "d",
    "maxPurity", "d",
    "minPloidy", "d",
    "maxPloidy", "d",
    "minDiploidProportion", "d",
    "maxDiploidProportion", "d",
    "version", "c",
    "somaticPenalty", "d",
    "wholeGenomeDuplication", "c",
    "msIndelsPerMb", "d",
    "msStatus", "c",
    "tml", "d",
    "tmlStatus", "c",
    "tmbPerMb", "d",
    "tmbStatus", "c",
    "svTumorMutationalBurden", "d",
    "runMode", "c",
    "targeted", "c"
  )

  ctypes <- paste(tab$type, collapse = "")
  purple_purity <- readr::read_tsv(x, col_types = ctypes)
  assertthat::assert_that(ncol(purple_purity) == nrow(tab))
  assertthat::assert_that(all(colnames(purple_purity) == tab$column))

  purple_purity <- purple_purity |>
    dplyr::mutate(
      dplyr::across(tidyselect::vars_select_helpers$where(is.numeric), round, 2),
      dplyr::across(dplyr::everything(), as.character)
    ) |>
    tidyr::pivot_longer(dplyr::everything(), names_to = "column", values_to = "value") |>
    dplyr::left_join(tab, by = "column") |>
    dplyr::mutate(value = toupper(.data$value)) |>
    dplyr::select("column", "value")

  p <- structure(purple_purity$value, names = purple_purity$column)

  summary <- dplyr::tribble(
    ~n, ~variable, ~value, ~details,
    2, "Purity", glue::glue('{p["purity"]} ({p["minPurity"]}-{p["maxPurity"]})'),
    "Purity of tumor in the sample (and min-max with score within 10% of best).",
    3, "Ploidy", glue::glue('{p["ploidy"]} ({p["minPloidy"]}-{p["maxPloidy"]})'),
    "Average ploidy of tumor sample after adjusting for purity (and min-max with score within 10% of best).",
    4, "Gender", glue::glue('{p["gender"]}'),
    "Gender as inferred by AMBER/COBALT.",
    7, "WGD", glue::glue('{p["wholeGenomeDuplication"]}'),
    "Whole genome duplication (more than 10 autosomes have average major allele ploidy > 1.5).",
    8, "MSI (indels/Mb)", glue::glue('{p["msStatus"]} ({p["msIndelsPerMb"]})'),
    "MSI status (MSI, MSS or UNKNOWN if somatic variants not supplied) & MS Indels per Mb.",
    9, "PolyclonalProp", glue::glue('{p["polyclonalProportion"]}'),
    "Proportion of CN regions that are more than 0.25 from a whole CN",
    10, "DiploidyProp", glue::glue('{p["diploidProportion"]} ({p["minDiploidProportion"]}-{p["maxDiploidProportion"]})'),
    "Proportion of CN regions that have 1 (+- 0.2) minor and major allele.",
    11, "TMB", glue::glue('{p["tmbPerMb"]} ({p["tmbStatus"]})'),
    paste(
      "Tumor mutational burden (# PASS variants per Megabase)",
      "(Status: 'HIGH' (>10 PASS per Mb), 'LOW' or 'UNKNOWN')."
    ),
    12, "TML", glue::glue('{p["tml"]} ({p["tmlStatus"]})'),
    "Tumor mutational load (# of missense variants) (Status: 'HIGH', 'LOW' or 'UNKNOWN').",
    13, "TMB-SV", glue::glue('{p["svTumorMutationalBurden"]}'),
    "# of non inferred, non single passing SVs."
  )

  list(
    raw = purple_purity,
    summary = summary
  )
}

#' Read PURPLE Somatic SNV VCF
#'
#' Reads the `purple.somatic.vcf.gz` file.
#'
#' @param x Path to the `purple.somatic.vcf.gz` file.
#'
#' @return A list with the input file as a tibble (with specific INFO fields)
#' and a tibble with a description of each INFO field.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.somatic.vcf.gz", package = "gpgr")
#' (snv <- purple_snv_vcf_read(x))
#' @export
purple_snv_vcf_read <- function(x) {
  assertthat::assert_that(file.exists(x), is_vcf(x))
  d <- bedr::read.vcf(x, split.info = TRUE, verbose = FALSE)
  info <-
    tibble::as_tibble(d$header[["INFO"]]) |>
    dplyr::select("ID", "Description")

  info_cols <- c(
    "PURPLE_AF", "PURPLE_CN",
    "PURPLE_GERMLINE", "PURPLE_MACN", "PURPLE_VCN",
    "HMF_HOTSPOT", "KT", "MH", "SUBCL", "TNC"
  )
  description <- dplyr::filter(info, .data$ID %in% info_cols)

  d <- tibble::as_tibble(d$vcf[c("CHROM", "POS", info_cols)])
  list(
    data = d,
    description = description
  )
}

#' Get PURPLE Kataegis Regions
#'
#' Reads the `purple.somatic.vcf.gz` file and extracts variants
#' within kataegis regions.
#'
#' @param x Path to the `purple.somatic.vcf.gz` file.
#'
#' @return A list with a tibble containing variants in kataegis clusters and
#' various metrics for each variant, and a tibble with a description of each
#' metric.
#'
#' @examples
#' x <- system.file("extdata/purple/purple.somatic.vcf.gz", package = "gpgr")
#' (k <- purple_kataegis(x))
#' @export
purple_kataegis <- function(x) {
  d <- purple_snv_vcf_read(x)
  info_cols <- c(
    "KT", "PURPLE_AF", "PURPLE_CN",
    "PURPLE_MACN", "PURPLE_VCN", "SUBCL",
    "MH", "TNC"
  )

  data <- d$data |>
    dplyr::filter(!is.na(.data$KT)) |>
    dplyr::select(c("CHROM", "POS", info_cols))

  description <- d$description |>
    dplyr::filter(.data$ID %in% info_cols) |>
    dplyr::arrange(.data$ID)

  list(
    data = data,
    description = description
  )
}
