#' Split Two-Field Column
#'
#' Splits a column with 2 comma-separated fields into
#' two columns.
#'
#' @param .data Input tibble data.
#' @param col Column to split.
#' @param is_pct Multiply by 100 (logical).
#'
#' @return A modified tibble with two columns.
#'
#' @examples
#' x <- tibble::tibble(a = letters[1:10], b = paste0(round(runif(10), 2), ",", round(runif(10), 2)))
#' (s <- gpgr:::split_double_col(x, b))
#'
#' @testexamples
#' expect_equal(colnames(s), c("a", "b1", "b2", "b"))
#' expect_error(gpgr:::split_double_col(x, c))
#' expect_equal(unname(sapply(s, class)), c("character", "numeric", "numeric", "numeric"))
#'
split_double_col <- function(.data, col, is_pct = FALSE) {
  # - separate field into two parts
  # - mutate to pct accordingly
  # - original field is mean of two parts
  f_q <- rlang::enquo(col)
  f_str <- rlang::quo_name(f_q)
  f1_str <- paste0(f_str, '1')
  f2_str <- paste0(f_str, '2')
  f1_q <- rlang::sym(f1_str)
  f2_q <- rlang::sym(f2_str)
  .data %>%
    tidyr::separate(!!f_q, c(f1_str, f2_str), sep = ",", fill = "right") %>%
    dplyr::mutate(
      !!f1_q := round(as.double(!!f1_q) * ifelse(is_pct, 100, 1), 1),
      !!f2_q := round(as.double(!!f2_q) * ifelse(is_pct, 100, 1), 1),
      !!f_q  := round(((!!f1_q + ifelse(is.na(!!f2_q), !!f1_q, !!f2_q)) / 2), 1)
    )
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

effect_abbreviations <- c(
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
effect_abbrev_nms <- names(effect_abbreviations)

abbreviate_effect <- function(effects) {
  # take string as x&y&z
  # split by &
  # abbreviate each piece and glue back with comma

  .abbreviate_effect <- function(effect) {
    ifelse(effect %in% effect_abbrev_nms, effect_abbreviations[effect], effect)
  }

  strsplit(effects, "&")[[1]] %>%
    purrr::map_chr(.abbreviate_effect) %>%
    paste(collapse = ", ")
}
abbreviate_effectv <- Vectorize(abbreviate_effect)

sv_path <- params$somatic_sv
sv_unmelted <- NULL
sv_all <- NULL

if (length(readLines(con = sv_path, n = 2)) > 1) {
  somatic_sv_tsv <- readr::read_tsv(sv_path, col_names = TRUE, col_types = "ccciicccccicccccdciicc")
  sv_unmelted <- somatic_sv_tsv %>%
    dplyr::select(-caller, -sample) %>%
    split_double_col(AF_BPI) %>%
    split_double_col(AF_PURPLE) %>%
    split_double_col(CN_PURPLE) %>%
    split_double_col(CN_change_PURPLE) %>%
    dplyr::mutate(
      AF_BPI = ifelse(is.na(AF_BPI), NA, paste0(AF_BPI, " (", AF_BPI1, ",", AF_BPI2, ")")),
      AF_PURPLE = ifelse(is.na(AF_PURPLE), NA, paste0(AF_PURPLE, " (", AF_PURPLE1, ", ", AF_PURPLE2, ")")),
      CN_PURPLE = ifelse(is.na(CN_PURPLE), NA, paste0(CN_PURPLE, " (", CN_PURPLE1, ", ", CN_PURPLE2, ")")),
      CN_change_PURPLE = ifelse(is.na(CN_change_PURPLE), NA, paste0(CN_change_PURPLE, " (", CN_change_PURPLE1, ", ", CN_change_PURPLE2, ")"))
    ) %>%
    tidyr::separate(split_read_support, c("SR_ref", "SR_alt"), ",", convert = TRUE) %>%
    tidyr::separate(paired_support_PR, c("PR_ref", "PR_alt"), ",", convert = TRUE) %>%
    dplyr::mutate(
      SR_PR_alt = paste0(SR_alt, ",", PR_alt),
      SR_PR_ref = paste0(SR_ref, ",", PR_ref),
      Ploidy = round(as.double(Ploidy_PURPLE), 2),
      chrom = sub("chr", "", chrom),
      chrom = as.factor(chrom), # for better DT filtering
      svtype = ifelse(is.na(PURPLE_status), svtype, "PURPLE_inf"),
      Start = ifelse(is.na(PURPLE_status), START_BPI, start),
      nann = count_pieces(annotation, ","),
      varnum = dplyr::row_number(),
      varnum = sprintf(glue::glue("%0{nchar(nrow(.))}d"), varnum))

  # Fix BND IDs
  sv_unmelted_bnd <- sv_unmelted %>%
    dplyr::filter(svtype == "BND") %>%
    tidyr::separate(ID, into = c("BND_group", "BND_mate"), sep = -1, convert = TRUE, remove = FALSE) %>%
    dplyr::group_by(BND_group) %>%
    dplyr::mutate(
      BND_ID = dplyr::group_indices(),
      BND_ID = sprintf(glue::glue("%0{nchar(nrow(.))}d"), BND_ID),
      BND_mate = ifelse(BND_mate == 0, "A", "B")) %>%
    dplyr::ungroup() %>%
    dplyr::bind_cols(., .[match(.$ID, .$MATEID), "chrom"]) %>%
    dplyr::rename(BND_mate_chrom = chrom1)

  sv_unmelted_other <- sv_unmelted %>%
    dplyr::filter(svtype != "BND")

  sv_unmelted <-
    dplyr::bind_rows(sv_unmelted_bnd,
                     sv_unmelted_other) %>%
    dplyr::mutate(
      END_BPI = base::format(END_BPI, big.mark = ",", trim = TRUE),
      End = ifelse(svtype == "BND", paste0(BND_mate_chrom, ":", END_BPI), END_BPI)) %>%
    dplyr::select(varnum, TierTop = tier,
                  Chr = chrom, Start, End,
                  Type = svtype,
                  ID, MATEID, BND_ID, BND_mate,
                  SR_PR_alt, SR_PR_ref, Ploidy,
                  AF_PURPLE, AF_BPI,
                  CNC = CN_change_PURPLE, CN = CN_PURPLE,
                  SScore = somaticscore, nann, annotation)

  sv_all <- sv_unmelted %>%
    dplyr::mutate(annotation = strsplit(annotation, ',')) %>%
    tidyr::unnest(annotation) %>%
    tidyr::separate(
      annotation, c('Event', 'Effect', 'Genes', 'Transcript', 'Detail', 'Tier'),
      sep = '\\|', convert = FALSE) %>%
    dplyr::mutate(
      ntrx = count_pieces(Transcript, "&"),
      ngen = count_pieces(Genes, "&"),
      neff = count_pieces(Effect, "&"),
      Transcript = Transcript %>% stringr::str_replace_all('&', ', '),
      Genes = Genes %>% stringr::str_replace_all('&', ', '),
      Effect = abbreviate_effectv(Effect)) %>%
    dplyr::distinct() %>%
    dplyr::arrange(Tier, Genes, Effect)
} else {
  sv_unmelted <- tibble(WARNING = "THERE WERE 0 SVs PRIORITISED!!")
  sv_all <- tibble(WARNING = "THERE WERE 0 SVs PRIORITISED!!")
}
