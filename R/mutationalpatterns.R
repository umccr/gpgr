#' Get COSMIC 2015 Signatures
#'
#' @usage data(cosmic_signatures_2015)
#' @docType data
#'
#' @format A matrix with 96 rows (one for each somatic mutation type within a specific context)
#' and 30 columns (one for each signature).
#'
"cosmic_signatures_2015"


#' Calculate COSMIC Signature Contribution
#'
#' Finds the linear combination of COSMIC (2015 and 2020) mutation signatures that
#' most closely reconstructs the SNV mutation matrix by solving the
#' nonnegative least-squares constraints problem.
#'
#' @param mut_mat Mutation count matrix (dimensions: m rows (mutation types) X 1 column (sample)).
#' @param signatures Signature matrix (dimensions: m rows (mutation types) X n columns (signatures))
#'
#' @return A list with the COSMIC 2015 and 2020 signature contributions to the
#' sample's signature.
#'
#' @export
sig_contribution <- function(mut_mat, signatures) {
  # Fit mutation matrix to cancer signatures
  fit_res <-
    MutationalPatterns::fit_to_signatures(mut_mat, signatures)$contribution %>%
    tibble::as_tibble(rownames = "sig") %>%
    dplyr::rename(contr = 2) %>%
    dplyr::filter(.data$contr > 0)

  if (nrow(fit_res) == 0) {
    fit_res1 <- tibble::tribble(
      ~sig, ~contr,
      "No Signatures found!", 0)
  }

  fit_res_contr <- fit_res %>%
    dplyr::mutate(
      contr = round(.data$contr, 0),
      RelFreq = round(.data$contr / sum(.data$contr), 2),
      Rank = as.integer(base::rank(-.data$contr))) %>%
    dplyr::select(.data$Rank, Signature = .data$sig,
                  Contribution = .data$contr, .data$RelFreq) %>%
    dplyr::arrange(.data$Rank)

  fit_res_contr
}


#' Count SNV Contexts
#'
#' Counts SNV Contexts.
#'
#' @param vcf_gr GRanges containing all mutation types from a single sample.
#' @param ref_genome The BSGenome reference genome to use.
#'
#' @return A list with two elements:
#' - snv_counts: matrix containing the number of SNVs per COSMIC context per gr.
#' - gr_snv: GRanges object containing the SNVs.
#'
#' @export
sig_count_snv <- function(vcf_gr, ref_genome) {
  gr_snv <- MutationalPatterns::get_mut_type(vcf_list = vcf_gr, type = "snv")
  snv_counts <- MutationalPatterns::mut_matrix(vcf_list = gr_snv, ref_genome = ref_genome)
  list(
    snv_counts = snv_counts,
    gr_snv = gr_snv
  )
}

#' Plot SNV Mutation Characteristics
#'
#' Plots SNV mutation characteristics.
#'
#' @param gr_snv GRanges containing SNVs from a single sample.
#' @param snv_counts A matrix with counts of SNV contexts.
#' @param ref_genome The BSGenome reference genome to use.
#'
#' @return A list with four ggplot2 objects:
#' - p_heatmap: a SNV mutation matrix as a heatmap.
#'   This is especially usefull when looking at a wide mutational context.
#' - p_river: a SNV mutation matrix as a riverplot.
#'   This is especially usefull when looking at a wide mutational context.
#' - p_96_profile: relative contribution of 96 trinucleotides.
#' - p_spectrum: point mutation spectrum.
#'
#' @export
sig_plot_snv <- function(gr_snv, snv_counts, ref_genome) {
  mut_to <- MutationalPatterns::mut_type_occurrences(
    vcf_list = gr_snv, ref_genome = ref_genome)

  mut_mat_ext_context <- MutationalPatterns::mut_matrix(
    vcf_list = gr_snv, ref_genome = ref_genome, extension = 2)

  p_spectrum <- MutationalPatterns::plot_spectrum(type_occurrences = mut_to, CT = TRUE,
                                                  condensed = TRUE, error_bars = "none") +
    ggplot2::theme(legend.position = "top")
  p_96_profile <- MutationalPatterns::plot_96_profile(mut_matrix = snv_counts, condensed = TRUE)
  p_heatmap <- MutationalPatterns::plot_profile_heatmap(mut_matrix = mut_mat_ext_context) +
    ggplot2::theme(legend.position = "none")
  p_river <- MutationalPatterns::plot_river(mut_matrix = mut_mat_ext_context) +
    ggplot2::theme(legend.position = "none")

  list(
    p_heatmap = p_heatmap,
    p_river = p_river,
    p_96_profile = p_96_profile,
    p_spectrum = p_spectrum
  )
}

#' Count INDEL Contexts
#'
#' Counts INDEL Contexts.
#'
#' @param vcf_gr GRanges containing all mutation types from a single sample.
#' @param ref_genome The BSGenome reference genome to use.
#'
#' @return A tibble containing the number of INDELs per COSMIC context per gr.
#'
#' @export
sig_count_indel <- function(vcf_gr, ref_genome) {
  gr_indel <- MutationalPatterns::get_mut_type(vcf_list = vcf_gr, type = "indel")
  gr_indel <- MutationalPatterns::get_indel_context(vcf_list = gr_indel, ref_genome = ref_genome)
  indel_counts <- MutationalPatterns::count_indel_contexts(vcf_list = gr_indel)
  indel_counts
}


#' Plot INDEL Mutation Characteristics
#'
#' Plots INDEL mutation characteristics.
#'
#' @param indel_counts INDEL context counts.
#'
#' @return A list with two ggplot2 objects:
#' - p_indel_main: the main INDEL contexts.
#' - p_indel_cont: the INDEL contexts.
#'
#' @export
sig_plot_indel <- function(indel_counts) {

  p_indel_main <- MutationalPatterns::plot_main_indel_contexts(counts = indel_counts) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.9, hjust = 1)) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank())
  p_indel_cont <- MutationalPatterns::plot_indel_contexts(counts = indel_counts, condensed = TRUE) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, vjust = 0.9, hjust = 1),
                   legend.position="top")

  list(
    p_indel_main = p_indel_main,
    p_indel_cont = p_indel_cont
  )
}

#' Count DBS Contexts
#'
#' Counts DBS Contexts.
#'
#' @param vcf_gr GRanges containing all mutation types from a single sample.
#'
#' @return A tibble containing the number of DBS per COSMIC context per gr.
#'
#' @export
sig_count_dbs <- function(vcf_gr) {
  gr_dbs <- MutationalPatterns::get_mut_type(vcf_list = vcf_gr, type = "dbs")
  gr_dbs <- MutationalPatterns::get_dbs_context(vcf_list = gr_dbs)
  dbs_counts <- MutationalPatterns::count_dbs_contexts(vcf_list = gr_dbs)
  dbs_counts
}

#' Plot DBS Mutation Characteristics
#'
#' Plots DBS mutation characteristics.
#'
#' @param dbs_counts DBS context counts.
#'
#' @return A list with two ggplot2 objects:
#' - p_dbs_main: the main DBS contexts.
#' - p_dbs_cont: the DBS contexts.
#'
#' @export
sig_plot_dbs <- function(dbs_counts) {

  p_dbs_main <- MutationalPatterns::plot_main_dbs_contexts(counts = dbs_counts)
  p_dbs_cont <- MutationalPatterns::plot_dbs_contexts(counts = dbs_counts, condensed = TRUE)

  list(
    p_dbs_main = p_dbs_main,
    p_dbs_cont = p_dbs_cont
  )
}

sig_contribution_table <- function() {
  sig_table <-
    readr::read_tsv(file = "misc/sig/v2_mar2015/signatures_description.tsv", col_types = "cc") %>%
    dplyr::mutate(Plot = paste0("![](misc/sig/v2_mar2015/img/sig-", signature, ".png)"),
                  signature = paste0("Sig", signature)) %>%
    dplyr::select(Signature = signature, Description = description, Plot)

  mut_sig_contr1 %>%
    dplyr::left_join(sig_table, by = "Signature") %>%
    knitr::kable() %>%
    kableExtra::kable_styling(c("hover", "striped"), font_size = 12) %>%
    kableExtra::scroll_box(height = "400px")

  sig_table2 <-
    readr::read_tsv(file = "misc/sig/v3_may2019/signatures_description.tsv", col_types = "cc") %>%
    dplyr::mutate(Plot = paste0("![](misc/sig/v3_may2019/img/sbs", signature, ".png)"),
                  signature = paste0("SBS", signature)) %>%
    dplyr::select(Signature = signature, Description = description, Plot)

  possible_seq_artefacts <- c("SBS27", "SBS43", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", "SBS51", "SBS52",
                              "SBS53", "SBS54", "SBS55", "SBS56", "SBS57", "SBS58", "SBS59", "SBS60")

  mut_sig_contr2 %>%
    dplyr::left_join(sig_table2, by = "Signature") %>%
    dplyr::mutate(Signature = ifelse(Signature %in% possible_seq_artefacts, paste0(Signature, " (SA)"), Signature)) %>%
    knitr::kable() %>%
    kableExtra::kable_styling(c("hover", "striped"), font_size = 12) %>%
    kableExtra::scroll_box(height = "400px")
}
