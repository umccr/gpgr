#' Split Two-Field Column
#'
#' Splits a column with 2 comma-separated fields into
#' two columns.
#'
#' @param d Input tibble data.
#' @param nms A character vector of columns to split.
#'
#' @return A modified tibble with the original columns averaged and split nicely.
#'
#' @examples
#' x <- tibble::tibble(a = letters[1:11],
#'                     b = c("0.4,0.8", paste0(round(runif(10), 2), ",", round(runif(10), 2))),
#'                     nacol = rep(NA, 11),
#'                     namix = sample(c(NA, "0.4,0.6"), 11, replace = TRUE))
#' (b <- gpgr:::split_double_col(x, "b"))
#' (nacol <- gpgr:::split_double_col(x, "nacol"))
#' (namix <- gpgr:::split_double_col(x, "namix"))
#'
#' @testexamples
#' expect_equal(colnames(b), "b")
#' expect_equal(nrow(x), 11)
#' expect_error(gpgr:::split_double_col(x, "c"))
#' expect_equal(b$b[1], "0.6 (0.4, 0.8)")
#'
split_double_col <- function(d, nms) {
  assertthat::assert_that(inherits(d, "tbl_df"))
  assertthat::assert_that(is.character(nms))

  outd <- d |>
    dplyr::select(nms) |>
    dplyr::mutate(num = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = nms, names_to = "col_nm", values_to = "x1_x2") |>
    tidyr::separate(.data$x1_x2, into = c("x1", "x2"), sep = ",", fill = "right") |>
    dplyr::mutate(x1 = round(as.double(.data$x1), 2),
                  x2 = round(as.double(.data$x2), 2))

  outd |>
    dplyr::mutate(avg = round(rowMeans(dplyr::select(outd, .data$x1, .data$x2), na.rm = TRUE), 2),
                  out = as.character(glue::glue("{avg} ({x1}, {x2})")),
                  out = ifelse(.data$out == "NaN (NA, NA)", NA_character_, .data$out)) |>
    dplyr::select(.data$num, .data$col_nm, .data$out) |>
    tidyr::pivot_wider(names_from = "col_nm", values_from = "out") |>
    dplyr::select(-.data$num)
}


#' Count Number of Parts in a String
#'
#' Counts number of pieces of a string separated by a pattern.
#' If it's an empty string, returns 0. If the pattern isn't found, returns 1.
#' If the pattern is found once, returns 2 (two pieces), etc.
#'
#' @param x Input string.
#' @param sep Pattern to count for.
#'
#' @return Number of parts.
#'
#' @examples
#' (a <- gpgr:::count_pieces("foo,bar,baz", sep = ","))
#' (b <- gpgr:::count_pieces("foo", sep = ","))
#' (k <- gpgr:::count_pieces("", sep = ","))
#' (m <- gpgr:::count_pieces(",", sep = ","))
#'
#'
#' @testexamples
#' expect_equal(a, 3)
#' expect_equal(b, 1)
#' expect_equal(k, 0)
#' expect_equal(m, 2)
#' expect_error(gpgr:::count_pieces("foo", NA))
#'
count_pieces <- function(x, sep) {
  ifelse(nchar(x) == 0, 0, stringr::str_count(x, sep) + 1)
}

#' Abbreviations for SV effects
#'
#' @usage data(EFFECT_ABBREVIATIONS)
#' @docType data
#'
#' @format Named character vector where names are the original names, and values
#' are the abbreviated names.
#'
"EFFECT_ABBREVIATIONS"

#' Abbreviate SV Effects
#'
#' Abbreviates SV effects column.
#'
#'
#' @param effects Input string with effects.
#'
#' @return An abbreviated string.
#'
#' @examples
#' (e1 <- gpgr:::abbreviate_effect("3_prime_UTR_truncation&start_lost&splice_region_variant"))
#' (e2 <- gpgr:::abbreviate_effect("duplication&foo&gene_fusion&BOOM&intron_variant"))
#' (e3 <- gpgr:::abbreviate_effect("TF_binding_site_variant&TFBS_ablation"))
#' (e4 <- gpgr:::abbreviate_effect("foo&bar&stop_gained&badaboom"))
#'
#' @testexamples
#' expect_equal(e1, "3UTRtrunc, SpliceRegV, StartLoss")
#' expect_equal(e2, "BOOM, Dup, foo, FusG, IntronV")
#' expect_equal(e3, "TFBSDel, TFBSVar")
#' expect_equal(e4, "badaboom, bar, foo, StopGain")
#'
abbreviate_effect <- function(effects) {

  effect_abbrev_nms <- names(gpgr::EFFECT_ABBREVIATIONS)

  # take string as x&y&z
  # split by &
  # abbreviate each piece and glue back with comma
  .abbreviate_effect <- function(effect) {
    ifelse(effect %in% effect_abbrev_nms, gpgr::EFFECT_ABBREVIATIONS[effect], effect)
  }

  strsplit(effects, "&")[[1]] |>
    purrr::map_chr(.abbreviate_effect) |>
    stringr::str_sort() |>
    paste(collapse = ", ")
}

#' Read SV TSV
#'
#' Reads the Manta TSV file output by umccrise.
#'
#' @param x Path to `manta.tsv` output by umccrise.
#'
#' @return A tibble corresponding to the input TSV file.
#'
#' @examples
#' x <- system.file("extdata/umccrise/sv/manta.tsv", package = "gpgr")
#' (sv <- umccrise_read_sv_tsv(x)$data)
#'
#' @testexamples
#' expect_equal(colnames(sv)[ncol(sv)], "ALT")
#'
#' @export
umccrise_read_sv_tsv <- function(x) {

  # tsv column names, description and types
  tab <- dplyr::tribble(
    ~Column, ~Description, ~Type,
    "caller", "Manta SV caller", "c",
    "sample", "Tumor sample name", "c",
    "chrom", "CHROM column in VCF", "c",
    "start", "POS column in VCF", "i",
    "end", "INFO/END: End position of the variant described in this record", "i",
    "svtype", "INFO/SVTYPE: Type of structural variant", "c",
    "split_read_support", "FORMAT/SR of tumor sample: Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999", "c",
    "paired_support_PE", "FORMAT/PE of tumor sample: ??", "c",
    "paired_support_PR", "FORMAT/PR of tumor sample: Spanning paired-read support for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999", "c",
    "AF_BPI", "INFO/BPI_AF: AF at each breakpoint (so AF_BPI1,AF_BPI2)", "c",
    "somaticscore", "INFO/SOMATICSCORE: Somatic variant quality score", "i",
    "tier", "INFO/SV_TOP_TIER (or 4 if missing): Highest priority tier for the effects of a variant entry", "c",
    "annotation", "INFO/SIMPLE_ANN: Simplified structural variant annotation: 'SVTYPE | EFFECT | GENE(s) | TRANSCRIPT | PRIORITY (1-4)'", "c",
    "AF_PURPLE", "INFO/PURPLE_AF: AF at each breakend (purity adjusted) (so AF_PURPLE1,AF_PURPLE2)", "c",
    "CN_PURPLE", "INFO/PURPLE_CN: CN at each breakend (purity adjusted) (so CN_PURPLE1,CN_PURPLE2)", "c",
    "CN_change_PURPLE", "INFO/PURPLE_CN_CHANGE: change in CN at each breakend (purity adjusted) (so CN_change_PURPLE1,CN_change_PURPLE2)", "c",
    "Ploidy_PURPLE", "INFO/PURPLE_PLOIDY: Ploidy of variant (purity adjusted)", "d",
    "PURPLE_status", "INFERRED if FILTER=INFERRED, or RECOVERED if has INFO/RECOVERED, else blank. INFERRED: Breakend inferred from copy number transition", "c",
    "START_BPI", "INFO/BPI_START: BPI adjusted breakend location", "i",
    "END_BPI", "INFO/BPI_END: BPI adjusted breakend location", "i",
    "ID", "ID column in VCF", "c",
    "MATEID", "INFO/MATEID: ID of mate breakend", "c",
    "ALT", "ALT column in VCF", "c")

  ctypes <- paste(tab$Type, collapse = "")
  somatic_sv_tsv <- readr::read_tsv(x, col_names = TRUE, col_types = ctypes)
  assertthat::assert_that(ncol(somatic_sv_tsv) == nrow(tab))
  assertthat::assert_that(all(colnames(somatic_sv_tsv) == tab$Column))
  list(
    data = somatic_sv_tsv,
    descr = tab
  )
}

#' Process Structural Variants
#'
#' Processes the Manta TSV file output by umccrise.
#'
#' @param x Path to `manta.tsv` output by umccrise.
#' @return A list with melted/unmelted tibbles (these are NULL if TSV file was empty).
#'
#' @examples
#' x <- system.file("extdata/umccrise/sv/manta.tsv", package = "gpgr")
#' (sv <- process_sv(x))
#'
#' @testexamples
#' expect_true(inherits(sv, "list"))
#' expect_equal(length(sv), 4)
#' expect_equal(names(sv), c("unmelted", "melted", "tsv_descr", "col_descr"))
#'
#' @export
process_sv <- function(x) {
  sv <- umccrise_read_sv_tsv(x)
  tsv_descr <- sv$descr
  sv <- sv$data

  if (nrow(sv) == 0) {
    return(list(
      unmelted = NULL,
      melted = NULL
    ))
  }
  col_descr <- dplyr::tribble(
    ~Column, ~Description,
    "nrow", "Row number that connects variants between tables in same tab set.",
    "vcfnum", "Original event row number that connects variants to events.",
    "TierTop", "Top priority of the event (from simple_sv_annotation: 1 highest, 4 lowest).",
    "Tier", "Priority of the specific event (from simple_sv_annotation: 1 highest, 4 lowest).",
    "Tier (Top)", "Priority of the specific event (Top tier of all event's annotations): 1 highest, 4 lowest. ",
    "Chr", "Chromosome.",
    "Start", "Start position as inferred by BPI. For PURPLE-inferred SVs we use POS.",
    "End", "End position. For BNDs = Chr:Start of the mate.Values are inferred by BPI (PURPLE-inferred SVs do not have an End).",
    "ID", "ID of BND from Manta. For PURPLE-inferred SVs this is PURPLE.",
    "MATEID", "ID of BND mate from Manta.",
    "BND_ID", "ID of BND pair simplified. BNDs with the same BND_ID belong to the same translocation event.",
    "BND_mate", "A or B depending on whether the first or second mate of the BND pair.",
    "Genes", "Genes involved in the event. DEL/DUP/INS events involving more than 2 genes are shown within separate table.",
    "Transcript", "Transcripts involved in the event. DEL/DUP/INS events involving more than 2 transcripts are shown within separate table.",
    "Effect", "SV effect (based on http://snpeff.sourceforge.net/SnpEff_manual.html#input).",
    "Detail", "Prioritisation detail (from simple_sv_annotation).",
    "Ploidy", "Ploidy of variant from PURPLE (purity adjusted).",
    "AF_PURPLE", "PURPLE AF at each breakend preceded by their average.",
    "AF_BPI", "BPI AF at each breakend preceded by their average.",
    "CN", "Copy Number at each breakend preceded by their average.",
    "CNC", "Copy Number Change at each breakend preceded by their average.",
    "SR_alt", "Number of Split Reads supporting the alt allele, where P(allele|read)>0.999.",
    "PR_alt", "Number of Paired Reads supporting the alt allele, where P(allele|read)>0.999.",
    "SR_PR_ref", "Number of Split Reads and Paired Reads supporting the ref allele, where P(allele|read)>0.999.",
    "SR_PR_sum", "Sum of SR_alt and PR_alt.",
    "Type", "Type of structural variant.",
    "SScore", "Somatic variant quality score.",
    "ntrx", "Number of transcripts for given event.",
    "ngen", "Number of genes for given event.",
    "nann", "Number of annotations for given event."
  )

  cols_to_split <- c("AF_BPI", "AF_PURPLE", "CN_PURPLE", "CN_change_PURPLE")
  double_cols <- split_double_col(sv, cols_to_split)
  unmelted <- sv |>
    dplyr::select(-dplyr::all_of(c(cols_to_split, "caller", "sample"))) |>
    dplyr::bind_cols(double_cols) |>
    tidyr::separate(.data$split_read_support, c("SR_ref", "SR_alt"), ",", convert = TRUE) |>
    tidyr::separate(.data$paired_support_PR, c("PR_ref", "PR_alt"), ",", convert = TRUE)
  unmelted <- unmelted |>
    dplyr::mutate(
      SR_PR_ref = paste0(.data$SR_ref, ",", .data$PR_ref),
      Ploidy = round(as.double(.data$Ploidy_PURPLE), 2),
      chrom = sub("chr", "", .data$chrom),
      svtype = ifelse(is.na(.data$PURPLE_status), .data$svtype, "PURPLE_inf"),
      Start = ifelse(is.na(.data$PURPLE_status), .data$START_BPI, .data$start),
      nann = count_pieces(.data$annotation, ","),
      vcfnum = dplyr::row_number(),
      vcfnum = sprintf(glue::glue("%0{nchar(nrow(unmelted))}d"), .data$vcfnum)) |>
    dplyr::rowwise() |>
    dplyr::mutate(SR_PR_sum = sum(.data$SR_alt, .data$PR_alt, na.rm = TRUE)) |>
    dplyr::ungroup()

  # BND IDs
  # Two BND mates share the same ID up to the last digit (0 or 1)
  unmelted_bnd1 <- unmelted |>
    dplyr::filter(.data$svtype == "BND") |>
    tidyr::separate(.data$ID, into = c("BND_group", "BND_mate"), sep = -1, convert = TRUE, remove = FALSE) |>
    dplyr::group_by(.data$BND_group)
  unmelted_bnd1 <- unmelted_bnd1 |>
    dplyr::mutate(
      # index per group 1, 2, 3..
      BND_ID = dplyr::cur_group_id(),
      # turns into 001, 002, 003... if you've got 100+ rows
      BND_ID = sprintf(glue::glue("%0{nchar(nrow(unmelted_bnd1))}d"), .data$BND_ID),
      BND_mate = ifelse(.data$BND_mate == 0, "A", "B")) |>
    dplyr::ungroup()

  # Grab each BND mate's chrom
  # Orphan mates have that info in the ALT field
  match_id2mateid <- match(unmelted_bnd1$ID, unmelted_bnd1$MATEID)
  unmelted_bnd2 <- unmelted_bnd1[match_id2mateid, c("chrom")] |>
    dplyr::rename(BND_mate_chrom = .data$chrom)

  unmelted_bnd <- dplyr::bind_cols(unmelted_bnd1, unmelted_bnd2) |>
    dplyr::mutate(
      BND_mate_chrom = ifelse(is.na(.data$BND_mate_chrom),
                              sub(".*chr(.*):.*", "orphan_\\1", .data$ALT),
                              .data$BND_mate_chrom))

  unmelted_other <- unmelted |>
    dplyr::filter(.data$svtype != "BND")

  unmelted_all <-
    dplyr::bind_rows(unmelted_bnd,
                     unmelted_other) |>
    dplyr::mutate(
      END_BPI = base::format(.data$END_BPI, big.mark = ",", trim = TRUE),
      Start = base::format(.data$Start, big.mark = ",", trim = TRUE),
      End = paste0(
        ifelse(.data$svtype == "BND", .data$BND_mate_chrom, .data$chrom),
        ":",
        .data$END_BPI),
      Start = paste0(.data$chrom, ":", .data$Start)) |>
    dplyr::select(.data$vcfnum, .data$nann, TierTop = .data$tier,
                  .data$Start, .data$End,
                  Type = .data$svtype,
                  .data$BND_ID, .data$BND_mate,
                  .data$SR_alt, .data$PR_alt, .data$SR_PR_sum, .data$SR_PR_ref, .data$Ploidy,
                  .data$AF_PURPLE, .data$AF_BPI,
                  CNC = .data$CN_change_PURPLE, CN = .data$CN_PURPLE,
                  SScore = .data$somaticscore, .data$annotation)

  abbreviate_effectv <- Vectorize(abbreviate_effect)

  melted <- unmelted_all |>
    dplyr::mutate(annotation = strsplit(.data$annotation, ',')) |>
    tidyr::unnest(.data$annotation) |>
    tidyr::separate(
      .data$annotation, c('Event', 'Effect', 'Genes', 'Transcript', 'Detail', 'Tier'),
      sep = '\\|', convert = FALSE) |>
    dplyr::mutate(
      ntrx = count_pieces(.data$Transcript, "&"),
      ngen = count_pieces(.data$Genes, "&"),
      neff = count_pieces(.data$Effect, "&"),
      Transcript = .data$Transcript |> stringr::str_replace_all('&', ', '),
      Genes = .data$Genes |> stringr::str_replace_all('&', ', '),
      Effect = abbreviate_effectv(.data$Effect),
      `Tier (Top)` = glue::glue("{Tier} ({TierTop})")) |>
    dplyr::distinct() |>
    dplyr::arrange(.data$`Tier (Top)`, .data$Genes, .data$Effect)

  list(
    unmelted = unmelted_all,
    melted = melted,
    tsv_descr = tsv_descr,
    col_descr = col_descr
  )
}

#' Line plot for SR, PR and SR + PR for BNDs
#'
#' Plots the number of split reads (`SR`), paired end reads (`PR`), and their
#' sum (`tot`) across all BNDs, sorted by `tot`.
#'
#' @param d A data.frame with an SR_PR_alt column.
#' @param title Main title of plot.
#' @param subtitle Subtitle of plot.
#'
#' @return A ggplot2 plot object.
#'
#' @examples
#' x <- system.file("extdata/umccrise/sv/manta.tsv", package = "gpgr")
#' d <- process_sv(x)$unmelted
#' plot_bnd_sr_pr_tot_lines(d)
#'
#' @export
plot_bnd_sr_pr_tot_lines <- function(d,
                                     title = "SR, PR and SR + PR line plot for BNDs",
                                     subtitle = "Events are sorted by decreasing tot values.") {
  assertthat::assert_that(all(c("Type", "SR_alt", "PR_alt") %in% colnames(d)))
  dplot <- d |>
    dplyr::filter(.data$Type == "BND") |>
    dplyr::select(SR = .data$SR_alt, PR = .data$PR_alt, Tier = .data$TierTop, .data$BND_ID) |>
    dplyr::distinct() |>
    dplyr::mutate(PR = ifelse(is.na(.data$PR), 0, .data$PR),
                  SR = ifelse(is.na(.data$SR), 0, .data$SR)) |>
    dplyr::rowwise() |>
    dplyr::mutate(tot = sum(.data$SR, .data$PR, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(.data$tot)) |>
    dplyr::mutate(bnd_event = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = c(.data$SR, .data$PR, .data$tot),
                        names_to = "Metric", values_to = "Count")

  p_all <- dplot |>
    ggplot2::ggplot(ggplot2::aes(x = .data$bnd_event, y = .data$Count, colour = .data$Metric)) +
    ggplot2::geom_line(alpha = 0.5) +
    ggplot2::geom_point(alpha = 0.5) +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor.x = ggplot2::element_blank()) +
    ggplot2::labs(title = title, subtitle = subtitle)

  # handle cases where no BNDs were detected
  p_tier <- ggplot2::ggplot() +
    ggplot2::theme_void() +
    ggplot2::geom_text(ggplot2::aes(x = 0, y = 0, label = "No BNDs detected!")) +
    ggplot2::xlab(NULL)

  if (nrow(dplot) > 0) {
    p_tier <- p_all +
      ggplot2::facet_wrap(~Tier) +
      ggplot2::labs(title = NULL, subtitle = paste(subtitle, "Faceted by Tier."))
  }
  list(p_all = p_all,
       p_tier = p_tier)
}

#' Histogram for SR, PR and SR + PR for BNDs
#'
#' Plots histograms for the number of split reads (`SR`), paired end reads (`PR`), and their
#' sum (`tot`) across all BNDs. Observations where the SR or PR value is 0 (NA) are not shown.
#'
#' @param d A data.frame with an SR_PR_alt column.
#' @param title Main title of plot.
#' @param subtitle Subtitle of plot.
#'
#' @return A ggplot2 plot object.
#'
#' @examples
#' x <- system.file("extdata/umccrise/sv/manta.tsv", package = "gpgr")
#' d <- process_sv(x)$unmelted
#' plot_bnd_sr_pr_tot_hist(d, "a title")
#'
#' @export
plot_bnd_sr_pr_tot_hist <- function(d,
                                    title = "SR, PR and SR + PR histogram for BNDs",
                                    subtitle = "Values of 0 (NA) are not shown.") {
  assertthat::assert_that(all(c("Type", "SR_alt", "PR_alt") %in% colnames(d)))
  dplot <- d |>
    dplyr::filter(.data$Type == "BND") |>
    dplyr::select(SR = .data$SR_alt, PR = .data$PR_alt, .data$BND_ID) |>
    dplyr::distinct() |>
    dplyr::mutate(PR = ifelse(is.na(.data$PR), 0, .data$PR),
                  SR = ifelse(is.na(.data$SR), 0, .data$SR)) |>
    dplyr::rowwise() |>
    dplyr::mutate(tot = sum(.data$SR, .data$PR, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(cols = c(.data$SR, .data$PR, .data$tot),
                        names_to = "Metric", values_to = "Value")

  if (nrow(dplot) > 0) {
    dplot |>
      dplyr::filter(.data$Value > 0) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$Value, fill = .data$Metric)) +
      ggplot2::geom_histogram(binwidth = 1) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = title, subtitle = subtitle) +
      ggplot2::facet_wrap(~.data$Metric, ncol = 1, scales = "free_y")
  } else {
    # handle cases where no BNDs were detected
    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(x = 0, y = 0, label = "No BNDs detected!")) +
      ggplot2::xlab(NULL)
  }

}
