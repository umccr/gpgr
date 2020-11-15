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
#'
#' @return A list with the COSMIC 2015 and 2020 signature contributions to the
#' sample's signature.
#'
#' @export
sig_contribution <- function(mut_mat) {
  sigs_2015 <- gpgr::cosmic_signatures_2015
  sigs_2020 <- get_known_signatures(muttype = "snv",
                                    source = "COSMIC",
                                    incl_poss_artifacts = TRUE)

  # Fit mutation matrix to cancer signatures
  fit_res1 <-
    MutationalPatterns::fit_to_signatures(mut_mat, sigs_2015)$contribution %>%
    tibble::as_tibble(rownames = "sig") %>%
    dplyr::rename(contr = 2) %>%
    dplyr::filter(.data$contr > 0)
  fit_res2 <-
    MutationalPatterns::fit_to_signatures(mut_mat, sigs_2020)$contribution %>%
    tibble::as_tibble(rownames = "sig") %>%
    dplyr::rename(contr = 2) %>%
    dplyr::filter(.data$contr > 0)

  if (nrow(fit_res1) == 0) {
    fit_res1 <- tibble::tribble(
      ~sig, ~contr,
      "No Signatures found!", 0)
  }

  if (nrow(fit_res2) == 0) {
    fit_res2 <- tibble::tribble(
      ~sig, ~contr,
      "No Signatures found!", 0)
  }

  fit_res_contr1 <- fit_res1 %>%
    dplyr::mutate(
      contr = round(.data$contr, 0),
      RelFreq = round(.data$contr / sum(.data$contr), 2),
      Rank = as.integer(base::rank(-.data$contr))) %>%
    dplyr::select(.data$Rank, Signature = .data$sig,
                  Contribution = .data$contr, .data$RelFreq) %>%
    dplyr::arrange(.data$Rank)

  fit_res_contr2 <- fit_res2 %>%
    dplyr::mutate(
      contr = round(.data$contr, 0),
      RelFreq = round(.data$contr / sum(.data$contr), 2),
      Rank = as.integer(base::rank(-.data$contr))) %>%
    dplyr::select(.data$Rank, Signature = .data$sig, Contribution = .data$contr, .data$RelFreq) %>%
    dplyr::arrange(.data$Rank)

  list(
    cosmic_2015_fit = fit_res_contr1,
    cosmic_2020_fit = fit_res_contr2
  )
}

#' Plot SNV Mutation Characteristics
#'
#' Plots SNV mutation characteristics.
#'
#' @param vcf_gr GRanges containing all mutation types from a single sample.
#' @param ref_genome The BSGenome reference genome object.
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
sig_plot_snv <- function(vcf_gr, ref_genome) {
  gr_snv <- MutationalPatterns::get_mut_type(
    vcf_list = gr, type = "snv"
  )
  mut_to <- MutationalPatterns::mut_type_occurrences(
    vcf_list = gr_snv, ref_genome = ref_genome
  )
  mut_mat <- MutationalPatterns::mut_matrix(
    vcf_list = gr_snv, ref_genome = ref_genome
  )
  mut_mat_ext_context <- MutationalPatterns::mut_matrix(
    vcf_list = gr_snv, ref_genome = ref_genome, extension = 2
  )

  p_spectrum <- MutationalPatterns::plot_spectrum(
    type_occurrences = mut_to, CT = TRUE, condensed = TRUE, error_bars = "none") +
    ggplot2::theme(legend.position = "bottom")
  p_96_profile <- MutationalPatterns::plot_96_profile(mut_matrix = mut_mat, condensed = TRUE)
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
