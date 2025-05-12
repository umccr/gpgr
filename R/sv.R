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
#' x <- tibble::tibble(
#'   a = letters[1:11],
#'   b = c("0.4,0.8", paste0(round(runif(10), 2), ",", round(runif(10), 2))),
#'   nacol = rep(NA, 11),
#'   namix = sample(c(NA, "0.4,0.6"), 11, replace = TRUE)
#' )
#' (b <- gpgr:::split_double_col(x, "b"))
#' (nacol <- gpgr:::split_double_col(x, "nacol"))
#' (namix <- gpgr:::split_double_col(x, "namix"))
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
    dplyr::select(dplyr::all_of(nms)) |>
    dplyr::mutate(num = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = dplyr::all_of(nms), names_to = "col_nm", values_to = "x1_x2") |>
    tidyr::separate(.data$x1_x2, into = c("x1", "x2"), sep = ",", fill = "right") |>
    dplyr::mutate(
      x1 = round(as.double(.data$x1), 2),
      x2 = round(as.double(.data$x2), 2)
    )

  outd |>
    dplyr::mutate(
      avg = round(rowMeans(dplyr::select(outd, "x1", "x2"), na.rm = TRUE), 2),
      out = as.character(glue::glue("{avg} ({x1}, {x2})")),
      out = ifelse(.data$out == "NaN (NA, NA)", NA_character_, .data$out)
    ) |>
    dplyr::select("num", "col_nm", "out") |>
    tidyr::pivot_wider(names_from = "col_nm", values_from = "out") |>
    dplyr::select(-"num")
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
#' @testexamples
#' expect_equal(a, 3)
#' expect_equal(b, 1)
#' expect_equal(k, 0)
#' expect_equal(m, 2)
#' expect_error(gpgr:::count_pieces("foo", NA))
#'
#' @export
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
#' @testexamples
#' expect_equal(e1, "3UTRtrunc, SpliceRegV, StartLoss")
#' expect_equal(e2, "BOOM, Dup, foo, FusG, IntronV")
#' expect_equal(e3, "TFBSDel, TFBSVar")
#' expect_equal(e4, "badaboom, bar, foo, StopGain")
#'
#' @export
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

sash_read_sv_tsv <- function(x) {
  tab <- dplyr::tribble(
    ~Column,         ~Description,                                                                                       ~Type,
    "chrom",         "CHROM column in VCF",                                                                               "c",
    "start",         "POS column in VCF",                                                                                 "i",
    "svtype",        "INFO/SVTYPE: Type of structural variant",                                                           "c",
    "VF_alt",        "Alternate allele variant fraction (VF) from esvee",                                                 "f",
    "DF_alt",        "Count of discordant fragments with reads either side BE from esvee",                                "f",
    "SF_alt",        "Alternate allele split-read support fraction (SF) from esvee",                                      "f",
    "SR_ref",        "FORMAT/REF of tumor sample",                                                                        "i",
    "PR_ref",        "FORMAT/REFPAIR of tumor sample",                                                                    "i",
    "QUAL",          "QUAL column in VCF",                                                                                "f",
    "tier",          "INFO/SV_TOP_TIER (or 4 if missing): Highest priority tier for the effects of a variant entry",      "c",
    "annotation",    "INFO/SIMPLE_ANN: Simplified structural variant annotation: 'SVTYPE | EFFECT | GENE(s) | TRANSCRIPT | PRIORITY (1-4)'", "c",
    "AF_PURPLE",     "INFO/PURPLE_AF: AF at each breakend (purity adjusted) (so AF_PURPLE1,AF_PURPLE2)",                  "c",
    "CN_PURPLE",     "INFO/PURPLE_CN: CN at each breakend (purity adjusted) (so CN_PURPLE1,CN_PURPLE2)",                  "c",
    "CN_change_PURPLE", "INFO/PURPLE_CN_CHANGE: change in CN at each breakend (purity adjusted) (so CN_change_PURPLE1,CN_change_PURPLE2)", "c",
    "PURPLE_status", "INFERRED if FILTER=INFERRED, or RECOVERED if has INFO/RECOVERED, else blank. INFERRED: Breakend inferred from copy number transition", "c",
    "ID",            "ID column in VCF",                                                                                  "c",
    "MATEID",        "INFO/MATEID: ID of mate breakend",                                                                  "c",
    "ALT",           "ALT column from VCF used in split_svs",                                                             "c",
    "DF_alt",        "Alternate allele strand bias (SB) from esvee",                                                      "f"
  )

  ctypes <- paste(tab$Type, collapse = "")
  d <- readr::read_tsv(x, col_names = TRUE, col_types = ctypes)
  assertthat::assert_that(ncol(d) == nrow(tab))
  assertthat::assert_that(all(colnames(d) == tab$Column))
  list(
    data = d,
    descr = tab
  )
}

split_svs <- function(x) {
  bps_types <- c("BND", "DEL", "DUP", "INS", "INV")

  x.grouped <- x |>
    dplyr::group_by(
      record_type = ifelse(.data$Type %in% bps_types, "bps", "other")
    )

  keys <- x.grouped |>
    dplyr::group_keys() |>
    dplyr::pull("record_type")

  x.split <- x.grouped |>
    dplyr::group_split(.keep = FALSE) |>
    purrr::set_names(keys)

  list(
    bps = purrr::pluck(x.split, "bps"),
    other = purrr::pluck(x.split, "other")
  )
}

join_breakpoint_entries <- function(x) {
  # Group by GRIDSS identifier (clipping trailing h/o [h: High, o: lOwer])
  bps <- x |>
    tidyr::separate("ID", into = c("BND_group", "BND_mate"), sep = -1, convert = TRUE, remove = FALSE) |>
    dplyr::group_by(.data$BND_group)

  # Set a sequential breakpoint identifier
  bps_groups <- bps |> dplyr::n_groups()
  bps |>
    dplyr::mutate(
      # Assign a unique ID based on current group
      BND_ID = sprintf(paste0("%0", nchar(bps_groups), "d"), dplyr::cur_group_id()),
      BND_mate = ifelse(.data$BND_mate == "o", "A", "B"),
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      end_position = sub("^.*:(\\d+).*$", "\\1", .data$ALT) |>
        as.numeric() |>
        base::format(big.mark = ",", trim = TRUE),
      end_chrom = sub("^.*chr(.*):.*$", "\\1", .data$ALT),
      end = paste0(.data$end_chrom, ":", .data$end_position),
    ) |>
    dplyr::select(-c("end_position", "end_chrom"))
}

remove_gene_fusion_dups <- function(.data, columns) {
  # Order elements of multi-entry effect values for reliable comparison
  v.groups <- c("frameshift_variant&gene_fusion", "gene_fusion")
  v.effects_ordered <- sapply(.data$Effect, function(s) {
    c1 <- stringr::str_split(s, "&") |> unlist()
    paste0(sort(c1), collapse = "&")
  })

  if (all(v.groups %in% v.effects_ordered)) {
    .data |> dplyr::filter(.data$Effect != "gene_fusion")
  } else {
    .data
  }
}

filter_and_split_annotations_sv <- function(x) {
  filter_conditions <- list(
    # Empty Gene field
    x$Genes == "",
    # Only Ensembl identifiers in the Gene field
    stringr::str_split(x$Genes, "[&-]|, ") |> purrr::map_lgl(\(x) stringr::str_starts(x, "ENSG") |> all()),
    # Fusions that do not involve genes
    x$Effect == "Fus",
    # PURPLE inferred SVs
    x$Type == "PURPLE_inf"
  )

  x.grouped <- x |>
    dplyr::group_by(
      filter = ifelse(purrr::reduce(filter_conditions, `|`), "filter", "retain")
    ) |>
    dplyr::select(-c("Type", "Top Tier"))

  keys <- x.grouped |>
    dplyr::group_keys() |>
    dplyr::pull("filter")

  x.split <- x.grouped |>
    dplyr::group_split(.keep = FALSE) |>
    purrr::set_names(keys)

  list(
    retained = purrr::pluck(x.split, "retain"),
    filtered = purrr::pluck(x.split, "filter")
  )
}

set_many_transcripts_sv <- function(x) {
  # Set many transcripts
  x.tmp <- x |>
    dplyr::rowwise() |>
    dplyr::mutate(
      `Transcript count` = stringr::str_split(.data$Transcripts, ", ") |> unlist() |> unique() |> length()
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      many_transcripts = ifelse(.data$`Transcript count` > 2, "many_transcripts", "few_transcripts")
    )

  # Build the many transcripts table
  mt <- x.tmp |>
    dplyr::filter(.data$many_transcripts == "many_transcripts") |>
    # Remove unneeded columns and rename others
    dplyr::select(
      "Record ID",
      "Annotation ID",
      "Tier (top)",
      "Start",
      "End",
      "Transcript count",
      "Transcripts",
    )

  # Clear transcript where appropriate
  x.ready <- x.tmp |>
    dplyr::mutate(
      Transcripts = ifelse(
        .data$many_transcripts == "few_transcripts" | is.na(.data$many_transcripts),
        .data$Transcripts,
        paste0("Many transcripts (", .data$`Transcript count`, ")")
      )
    ) |>
    dplyr::select(-c("many_transcripts", "Transcript count"))

  list(
    many_transcripts = mt,
    sv = x.ready
  )
}

#' Process SV TSV
#'
#' @param x Path to SV TSV.
#' @return List of many things.
#'
#' @export
process_sv <- function(x) {
  # Read input and set column information
  sv.input <- sash_read_sv_tsv(x)

  # Early return if no data to process
  if (nrow(sv.input$data) == 0) {
    return()
  }

  # Prepare input
  sv.ready <- sv.input$data |>
    dplyr::mutate(
      "annotation_count" = count_pieces(.data$annotation, ","),
      "Top Tier" = .data$tier,
      start = paste(.data$chrom, base::format(.data$start, big.mark = ",", trim = TRUE), sep = ":"),
      Type = ifelse(is.na(.data$PURPLE_status), .data$svtype, "PURPLE_inf"),
      "Record ID" = dplyr::row_number()
    ) |>
    dplyr::select(-c(
      "chrom",
      "PURPLE_status",
      "tier",
      "svtype"
    ))

  # Split out breakpoints for merging
  sv.split <- split_svs(sv.ready)

  # Complete breakpoint records with corresponding mate, then combine with non-breakend records
  sv.bps <- join_breakpoint_entries(sv.split$bps)
  sv.tmp <- dplyr::bind_rows(sv.bps, sv.split$other)

  # Format some columns
  cols_to_split <- c("AF_PURPLE", "CN_PURPLE")
  double_cols <- split_double_col(sv.tmp, cols_to_split)
  sv.tmp <- sv.tmp |>
    dplyr::select(-c(dplyr::all_of(cols_to_split))) |>
    dplyr::bind_cols(double_cols)

  # Format a table for to be used as the SV Map
  sv.map <- sv.tmp |>
    dplyr::select(
      "Record ID",
      "Annotations" = "annotation_count",
      "Top Tier",
      "Start" = "start",
      "End" = "end",
      "Type",
      "Breakend ID" = "BND_ID",
      "Breakend Mate" = "BND_mate",
      "VF_alt",
      "DF_alt",
      "SF_alt",
      "PURPLE AF" = "AF_PURPLE",
      "PURPLE CN" = "CN_PURPLE"
    ) |>
    dplyr::arrange(.data$`Record ID`)

  # Melt annotations
  sv.melted_all <- sv.tmp |>
    # Split into individual annotations
    dplyr::mutate(annotation = strsplit(.data$annotation, ",")) |>
    # Convert annotation fields into columns
    tidyr::unnest("annotation") |>
    tidyr::separate_wider_delim(
      cols = "annotation", delim = "|",
      names = c("Event", "Effect", "Genes", "Transcripts", "Detail", "Tier")
    ) |>
    # Remove gene_fusion annotations for variants where frameshift_variant&gene_fusion already exist
    dplyr::group_by(dplyr::across(-"Effect")) |>
    dplyr::group_modify(remove_gene_fusion_dups) |>
    dplyr::ungroup() |>
    # Remove unused columns
    dplyr::select(-c("Event", "ALT")) |>
    # Create columns, modify others
    dplyr::mutate(
      "Annotation ID" = dplyr::row_number(),
      "Tier (top)" = paste0(.data$Tier, " (", .data$`Top Tier`, ")"),
      "Genes" = stringr::str_replace_all(.data$Genes, "&", ", "),
      "Transcripts" = stringr::str_replace_all(.data$Transcripts, "&", ", "),
    ) |>
    # Sort rows
    dplyr::arrange(.data$`Tier (top)`, .data$`Record ID`, .data$Genes, .data$Effect)

  # Abbreviate effects
  abbreviate_effectv <- Vectorize(abbreviate_effect)
  sv.melted_all$Effect <- abbreviate_effectv(sv.melted_all$Effect)

  # Select and rename columns
  sv.melted_all <- sv.melted_all |>
    dplyr::select(
      "Record ID",
      "Annotation ID",
      "Tier (top)",
      "Start" = "start",
      "End" = "end",
      "Breakend ID" = "BND_ID",
      "Breakend Mate" = "BND_mate",
      "Effect",
      "Genes",
      "Transcripts",
      "Effect",
      "Detail",
      "VF_alt",
      "DF_alt",
      "SF_alt",
      "PURPLE AF" = "AF_PURPLE",
      "PURPLE CN" = "CN_PURPLE",
      # Dropped after ops for non-map outputs
      "Top Tier",
      "Type"
    )

  # Create and set many transcript values
  sv.annotations.many_transcript_data <- set_many_transcripts_sv(sv.melted_all)
  sv.annotations.many_transcripts <- purrr::pluck(sv.annotations.many_transcript_data, "many_transcripts")
  sv.annotations <- purrr::pluck(sv.annotations.many_transcript_data, "sv")

  # Filter unwanted annotations
  sv.annotations.split <- filter_and_split_annotations_sv(sv.annotations)

  list(
    map = sv.map,
    map_melted = sv.melted_all,
    annotations = sv.annotations.split$retained,
    annotations_filtered = sv.annotations.split$filtered,
    many_transcripts = sv.annotations.many_transcripts
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
#' x <- system.file("extdata/sash/sv.prioritised.tsv", package = "gpgr")
#' d <- process_sv(x)$map
#' plot_bnd_sr_pr_tot_lines(d)
#' @export
plot_bnd_sr_pr_tot_lines <- function(d,
                                     title = "SR, PR and SR + PR line plot for BNDs",
                                     subtitle = "Events are sorted by decreasing tot values.") {
  assertthat::assert_that(all(c("Type", "SR_alt", "PR_alt") %in% colnames(d)))
  dplot <- d |>
    dplyr::filter(.data$Type == "BND") |>
    dplyr::select(SR = "SR_alt", PR = "PR_alt", Tier = "Top Tier", "Breakend ID") |>
    dplyr::distinct() |>
    dplyr::mutate(
      PR = ifelse(is.na(.data$PR), 0, .data$PR),
      SR = ifelse(is.na(.data$SR), 0, .data$SR)
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(tot = sum(.data$SR, .data$PR, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    dplyr::arrange(dplyr::desc(.data$tot)) |>
    dplyr::mutate(bnd_event = dplyr::row_number()) |>
    tidyr::pivot_longer(
      cols = c(.data$SR, .data$PR, .data$tot),
      names_to = "Metric", values_to = "Count"
    )

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
  list(
    p_all = p_all,
    p_tier = p_tier
  )
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
#' x <- system.file("extdata/sash/sv.prioritised.tsv", package = "gpgr")
#' d <- process_sv(x)$map
#' plot_bnd_sr_pr_tot_hist(d, "a title")
#' @export
plot_bnd_sr_pr_tot_hist <- function(d,
                                    title = "SR, PR and SR + PR histogram for BNDs",
                                    subtitle = "Values of 0 (NA) are not shown.") {
  assertthat::assert_that(all(c("Type", "SR_alt", "PR_alt") %in% colnames(d)))
  dplot <- d |>
    dplyr::filter(.data$Type == "BND") |>
    dplyr::select(SR = "SR_alt", PR = "PR_alt", "Breakend ID") |>
    dplyr::distinct() |>
    dplyr::mutate(
      PR = ifelse(is.na(.data$PR), 0, .data$PR),
      SR = ifelse(is.na(.data$SR), 0, .data$SR)
    ) |>
    dplyr::rowwise() |>
    dplyr::mutate(tot = sum(.data$SR, .data$PR, na.rm = TRUE)) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(
      cols = c(.data$SR, .data$PR, .data$tot),
      names_to = "Metric", values_to = "Value"
    )

  if (nrow(dplot) > 0) {
    dplot |>
      dplyr::filter(.data$Value > 0) |>
      ggplot2::ggplot(ggplot2::aes(x = .data$Value, fill = .data$Metric)) +
      ggplot2::geom_histogram(binwidth = 1) +
      ggplot2::theme_bw() +
      ggplot2::labs(title = title, subtitle = subtitle) +
      ggplot2::facet_wrap(~ .data$Metric, ncol = 1, scales = "free_y")
  } else {
    # handle cases where no BNDs were detected
    ggplot2::ggplot() +
      ggplot2::theme_void() +
      ggplot2::geom_text(ggplot2::aes(x = 0, y = 0, label = "No BNDs detected!")) +
      ggplot2::xlab(NULL)
  }
}
