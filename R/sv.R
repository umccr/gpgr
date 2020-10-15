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
#'                     namix = sample(c(NA, "0.4,0.6"), 11, replace = T))
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

  d %>%
    dplyr::select(nms) %>%
    dplyr::mutate(num = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = nms, names_to = "col_nm", values_to = "x1_x2") %>%
    tidyr::separate(.data$x1_x2, into = c("x1", "x2"), sep = ",", fill = "right") %>%
    dplyr::mutate(x1 = round(as.double(.data$x1), 2),
                  x2 = round(as.double(.data$x2), 2)) %>%
    dplyr::mutate(avg = round(rowMeans(dplyr::select(., .data$x1, .data$x2), na.rm = TRUE), 2),
                  out = as.character(glue::glue("{avg} ({x1}, {x2})")),
                  out = ifelse(out == "NaN (NA, NA)", NA_character_, out)) %>%
    dplyr::select(.data$num, .data$col_nm, .data$out) %>%
    tidyr::pivot_wider(names_from = "col_nm", values_from = "out") %>%
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

EFFECT_ABBREVIATIONS <- c(
  "3_prime_UTR_truncation" = "3UTRtrunc", "3_prime_UTR_variant" = "3UTRvar",
  "5_prime_UTR_truncation" = "5UTRtrunc",  "5_prime_UTR_variant" = "5UTRvar",
  "feature_fusion" = "Fus", "bidirectional_gene_fusion" = "BidFusG", "gene_fusion" = "FusG",
  "chromosome_number_variation" = "ChromNumV", "conservative_inframe_deletion" = "ConsInframeDel",
  "downstream_gene_variant" = "DnstreamGV", "upstream_gene_variant" = "UpstreamGV",
  "duplication" = "Dup", "exon_loss_variant" = "ExonLossV",
  "feature_ablation" = "DelG", "transcript_ablation" = "DelTx",
  "frameshift_variant" = "FrameshiftV", "intergenic_region" = "IntergenReg", "intragenic_variant" = "IntragenV",
  "intron_variant" = "IntronV",  "no_func_effect" = "NoFuncEff", "no_prio_effect" = "NoPrioEff",
  "non_coding_transcript_variant" = "NoncodTxV",
  "splice_acceptor_variant" = "SpliceAccV",
  "splice_donor_variant" = "SpliceDonV", "splice_region_variant" = "SpliceRegV",
  "start_lost" = "StartLoss", "stop_gained" = "StopGain", "stop_lost" = "StopLoss",
  "TF_binding_site_variant" = "TFBSVar",  "TFBS_ablation" = "TFBSDel")


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

  effect_abbrev_nms <- names(EFFECT_ABBREVIATIONS)

  # take string as x&y&z
  # split by &
  # abbreviate each piece and glue back with comma
  .abbreviate_effect <- function(effect) {
    ifelse(effect %in% effect_abbrev_nms, EFFECT_ABBREVIATIONS[effect], effect)
  }

  strsplit(effects, "&")[[1]] %>%
    purrr::map_chr(.abbreviate_effect) %>%
    stringr::str_sort() %>%
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
#' (sv <- umccrise_read_sv_tsv(x))
#'
#' @testexamples
#' expect_equal(colnames(sv)[ncol(sv)], "ALT")
#'
#' @export
umccrise_read_sv_tsv <- function(x) {

  # tsv column names + types
  nm <- c("caller" = "c", "sample" = "c",
          "chrom" = "c", "start" = "i", "end" = "i", "svtype" = "c",
          "split_read_support" = "c", "paired_support_PE" = "c", "paired_support_PR" = "c",
          "AF_BPI" = "c", "somaticscore" = "i", "tier" = "c", "annotation" = "c",
          "AF_PURPLE" = "c", "CN_PURPLE" = "c", "CN_change_PURPLE" = "c", "Ploidy_PURPLE" = "d",
          "PURPLE_status" = "c", "START_BPI" = "i", "END_BPI" = "i", "ID" = "c",
          "MATEID" = "c", "ALT" = "c")

  ctypes <- paste(nm, collapse = "")
  somatic_sv_tsv <- readr::read_tsv(x, col_names = TRUE, col_types = ctypes)
  assertthat::assert_that(ncol(somatic_sv_tsv) == length(nm))
  assertthat::assert_that(all(colnames(somatic_sv_tsv) == names(nm)))
  somatic_sv_tsv
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
#' expect_equal(length(sv), 2)
#' expect_equal(names(sv), c("unmelted", "melted"))
#'
#' @export
process_sv <- function(x) {
  sv <- umccrise_read_sv_tsv(x)

  if (nrow(sv) == 0) {
    return(list(
      unmelted = NULL,
      melted = NULL
    ))
  }

  cols_to_split <- c("AF_BPI", "AF_PURPLE", "CN_PURPLE", "CN_change_PURPLE")
  double_cols <- split_double_col(sv, cols_to_split)
  unmelted <- sv %>%
    dplyr::select(-dplyr::all_of(c(cols_to_split, "caller", "sample"))) %>%
    dplyr::bind_cols(double_cols) %>%
    tidyr::separate(.data$split_read_support, c("SR_ref", "SR_alt"), ",", convert = TRUE) %>%
    tidyr::separate(.data$paired_support_PR, c("PR_ref", "PR_alt"), ",", convert = TRUE) %>%
    dplyr::mutate(
      SR_PR_alt = paste0(.data$SR_alt, ",", .data$PR_alt),
      SR_PR_ref = paste0(.data$SR_ref, ",", .data$PR_ref),
      Ploidy = round(as.double(.data$Ploidy_PURPLE), 2),
      chrom = sub("chr", "", .data$chrom),
      svtype = ifelse(is.na(.data$PURPLE_status), .data$svtype, "PURPLE_inf"),
      Start = ifelse(is.na(.data$PURPLE_status), .data$START_BPI, .data$start),
      nann = count_pieces(.data$annotation, ","),
      vcfnum = dplyr::row_number(),
      vcfnum = sprintf(glue::glue("%0{nchar(nrow(.))}d"), .data$vcfnum))

  # BND IDs
  # Two BND mates share the same ID up to the last digit (0 or 1)
  unmelted_bnd1 <- unmelted %>%
    dplyr::filter(.data$svtype == "BND") %>%
    tidyr::separate(.data$ID, into = c("BND_group", "BND_mate"), sep = -1, convert = TRUE, remove = FALSE) %>%
    dplyr::group_by(.data$BND_group) %>%
    dplyr::mutate(
      # index per group 1, 2, 3..
      BND_ID = dplyr::cur_group_id(),
      # turns into 001, 002, 003... if you've got 100+ rows
      BND_ID = sprintf(glue::glue("%0{nchar(nrow(.))}d"), .data$BND_ID),
      BND_mate = ifelse(.data$BND_mate == 0, "A", "B")) %>%
    dplyr::ungroup()

  # Grab each BND mate's chrom
  # Orphan mates have that info in the ALT field
  match_id2mateid <- match(unmelted_bnd1$ID, unmelted_bnd1$MATEID)
  unmelted_bnd2 <- unmelted_bnd1[match_id2mateid, c("chrom")] %>%
    dplyr::rename(BND_mate_chrom = .data$chrom)

  unmelted_bnd <- dplyr::bind_cols(unmelted_bnd1, unmelted_bnd2) %>%
    dplyr::mutate(
      BND_mate_chrom = ifelse(is.na(BND_mate_chrom), sub(".*chr(.*):.*", "orphan_\\1", ALT), BND_mate_chrom))

  unmelted_other <- unmelted %>%
    dplyr::filter(.data$svtype != "BND")

  unmelted_all <-
    dplyr::bind_rows(unmelted_bnd,
                     unmelted_other) %>%
    dplyr::mutate(
      END_BPI = base::format(.data$END_BPI, big.mark = ",", trim = TRUE),
      Start = base::format(.data$Start, big.mark = ",", trim = TRUE),
      End = paste0(
        ifelse(.data$svtype == "BND", .data$BND_mate_chrom, .data$chrom),
        ":",
        .data$END_BPI),
      Start = paste0(.data$chrom, ":", .data$Start)) %>%
    dplyr::select(.data$vcfnum, TierTop = .data$tier,
                  .data$Start, .data$End,
                  Type = .data$svtype,
                  .data$ID, .data$MATEID, .data$BND_ID, .data$BND_mate,
                  .data$SR_PR_alt, .data$SR_PR_ref, .data$Ploidy,
                  .data$AF_PURPLE, .data$AF_BPI,
                  CNC = .data$CN_change_PURPLE, CN = .data$CN_PURPLE,
                  SScore = .data$somaticscore, .data$nann, .data$annotation)

  abbreviate_effectv <- Vectorize(abbreviate_effect)

  melted <- unmelted_all %>%
      dplyr::mutate(annotation = strsplit(.data$annotation, ',')) %>%
      tidyr::unnest(.data$annotation) %>%
      tidyr::separate(
        .data$annotation, c('Event', 'Effect', 'Genes', 'Transcript', 'Detail', 'Tier'),
        sep = '\\|', convert = FALSE) %>%
      dplyr::mutate(
        ntrx = count_pieces(.data$Transcript, "&"),
        ngen = count_pieces(.data$Genes, "&"),
        neff = count_pieces(.data$Effect, "&"),
        Transcript = .data$Transcript %>% stringr::str_replace_all('&', ', '),
        Genes = .data$Genes %>% stringr::str_replace_all('&', ', '),
        Effect = abbreviate_effectv(.data$Effect)) %>%
      dplyr::distinct() %>%
      dplyr::arrange(.data$Tier, .data$Genes, .data$Effect)

  list(
    unmelted = unmelted_all,
    melted = melted
  )
}

# histo + density plots of PR - SR for BNDs
plot_bnd_sr_pr <- function(d, nm) {
  assertthat::assert_that(all(c("Type", "SR_PR_alt") %in% colnames(d)))
  dplot <- d %>%
    dplyr::filter(Type == "BND") %>%
    dplyr::select(SR_PR_alt) %>%
    tidyr::separate(SR_PR_alt, into = c("SR", "PR"), convert = TRUE) %>%
    dplyr::mutate(PR = ifelse(is.na(PR), 0, PR),
                  SR = ifelse(is.na(SR), 0, SR),
                  PR_minus_SR = PR - SR)

  p1 <- dplot %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$PR_minus_SR)) +
    ggplot2::geom_histogram(fill = "darkblue", binwidth = 1) +
    ggplot2::theme_bw()
  p2 <- dplot %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$PR_minus_SR)) +
    ggplot2::geom_density(colour = "darkblue") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(nm)

  p1 / p2

}

# line plot for SR, PR and SR + PR for BNDs
plot_bnd_sr_pr_tot <- function(d, nm) {
  assertthat::assert_that(all(c("Type", "SR_PR_alt") %in% colnames(d)))
  dplot <- d %>%
    dplyr::filter(Type == "BND") %>%
    dplyr::select(SR_PR_alt) %>%
    tidyr::separate(SR_PR_alt, into = c("SR", "PR"), convert = TRUE) %>%
    dplyr::mutate(PR = ifelse(is.na(PR), 0, PR),
                  SR = ifelse(is.na(SR), 0, SR)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(tot = sum(SR, PR, na.rm = T)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(dplyr::desc(tot)) %>%
    dplyr::mutate(n = dplyr::row_number()) %>%
    tidyr::pivot_longer(cols = c(SR, PR, tot))

  p <- dplot %>%
    ggplot2::ggplot(ggplot2::aes(x = .data$n, y = .data$value, colour = .data$name)) +
    ggplot2::geom_line() +
    ggplot2::scale_x_continuous(breaks=scales::breaks_extended(10)) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(nm)

  p
}

# sv1 <- process_sv("~/Desktop/tmp/SBJ1_keep_pass.tsv")
# sv2 <- process_sv("~/Desktop/tmp/SBJ2_keep_pass.tsv")
# sv3 <- process_sv("~/Desktop/tmp/SEQC_keep_pass.tsv")
#
# p1 <- plot_bnd_sr_pr(sv1$unmelted, "574_1")
# p2 <- plot_bnd_sr_pr(sv2$unmelted, "574_2")
# p3 <- plot_bnd_sr_pr(sv3$unmelted, "SEQC")
# p1 | p2 | p3
#
# p1 <- plot_bnd_sr_pr_tot(sv1$unmelted, "574_1")
# p2 <- plot_bnd_sr_pr_tot(sv2$unmelted, "574_2")
# p3 <- plot_bnd_sr_pr_tot(sv3$unmelted, "SEQC")
# p1 | p2 | p3
