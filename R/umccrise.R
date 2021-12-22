#' HRDetect and CHORD Summary Table
#'
#' Return a summary table with HRDetect and CHORD results.
#'
#' @param chord_res Result from running [sigrap::chord_run()].
#' @param hrdetect_res Result from running [sigrap::hrdetect_run()].
#'
#' @return A list with a tibble and a gt_tbl object (see [gt::gt()]).
#'
#' @examples
#' \dontrun{
#' snv <- system.file("extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz", package = "gpgr")
#' sv <- system.file("extdata/umccrise/sv/manta.vcf.gz", package = "gpgr")
#' cnv <- system.file("extdata/purple/purple.cnv.somatic.tsv", package = "gpgr")
#' nm <- "SampleA"
#' genome <- "hg38"
#' snvoutdir <- tempdir()
#' hrdetect_res <- sigrap::hrdetect_run(nm, snv, sv, cnv, genome, snvoutdir)
#' chord_res <- sigrap::chord_run(
#'   vcf.snv = snv, sample.name = nm,
#'   df.sv = gpgr:::chord_mantavcf2df(sv)
#' )
#' hrd_results_tabs(hrdetect_res = hrdetect_res, chord_res = chord_res)
#' }
#'
#' @export
hrd_results_tabs <- function(hrdetect_res, chord_res) {
  sn <- chord_res$prediction[, "sample", drop = T]
  assertthat::are_equal(hrdetect_res[, "sample", drop = T], sn)

  hrdetect_res_tab <-
    tibble::tibble(
      col = colnames(hrdetect_res),
      val = unlist(hrdetect_res[1, ])
    ) |>
    dplyr::filter(.data$col != "sample")

  chord_res_tab <-
    tibble::tibble(
      col = colnames(chord_res$prediction),
      val = unlist(chord_res$prediction[1, ])
    ) |>
    dplyr::filter(col != "sample")


  colnames(hrdetect_res_tab) <- c("HRDetect", "results_hrdetect")
  colnames(chord_res_tab) <- c("CHORD", "results_chord")

  tab1 <-
    dplyr::bind_rows(
      hrdetect_res_tab,
      tibble::tribble(
        ~HRDetect, ~results_hrdetect,
        " ", " "
      )
    )
  tab2 <-
    dplyr::bind_rows(
      chord_res_tab[1:7, ],
      tibble::tribble(
        ~CHORD, ~results_chord,
        " ", " ",
        " ", " "
      )
    )
  tab3 <- chord_res_tab[8:16, ] |>
    purrr::set_names(c("CHORD2", "results_chord2"))

  hrd_results_tab <- dplyr::bind_cols(tab1, tab2, tab3)

  hrd_results_gt <-
    hrd_results_tab |>
    gt::gt() |>
    gt::tab_header(
      title = glue::glue("HRD Results for {sn}")
    ) |>
    gt::tab_spanner(
      label = "HRDetect",
      columns = c("HRDetect", "results_hrdetect")
    ) |>
    gt::tab_spanner(
      label = "CHORD",
      columns = c("CHORD", "results_chord", "CHORD2", "results_chord2")
    ) |>
    gt::cols_label(
      HRDetect = "",
      results_hrdetect = "",
      CHORD = "",
      results_chord = "",
      CHORD2 = "",
      results_chord2 = ""
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = c("HRDetect", "CHORD", "CHORD2")
      )
    ) |>
    gt::cols_align(
      align = "right",
      columns = c("results_hrdetect")
    ) |>
    gt::tab_style(
      style = gt::cell_borders(
        sides = "left",
        color = "#BBBBBB",
        weight = gt::px(2.5),
        style = "solid"
      ),
      locations = gt::cells_body(
        columns = c("CHORD"),
        rows = dplyr::everything()
      )
    ) |>
    gt::tab_options(table.align = "left")

  list(
    sample = sn,
    hrd_results_tab = hrd_results_tab,
    hrd_results_gt = hrd_results_gt
  )
}


#' umccrise AF summary
#'
#' Get summary table and plot for SNV allele frequencies output by umccrise.
#'
#' @param af_global_file Path to 'global' AF file.
#' @param af_keygenes_file Path to 'keygenes' AF file.
#'
#' @return A list containing a gt table and a plot summarising AFs.
#'
#' @export
af_summary <- function(af_global_file, af_keygenes_file) {
  af_global <-
    readr::read_tsv(af_global_file, col_types = "d") |>
    dplyr::mutate(set = "Global")

  af_keygenes <-
    readr::read_tsv(af_keygenes_file, col_types = "cicccd") |>
    dplyr::select(.data$af) |>
    dplyr::mutate(set = "Key genes CDS")

  af_both <-
    dplyr::bind_rows(af_global, af_keygenes) |>
    dplyr::mutate(set = factor(.data$set, levels = c("Global", "Key genes CDS")))

  mode2 <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }

  af_stats <- af_both |>
    dplyr::group_by(.data$set) |>
    dplyr::summarise(
      n = dplyr::n(),
      mean = round(base::mean(.data$af), 2),
      median = round(stats::median(.data$af), 2),
      mode = round(mode2(.data$af), 2),
      .groups = "drop_last"
    ) |>
    tidyr::complete(.data$set, fill = list(n = 0))

  af_stats_gt <- af_stats |>
    gt::gt(rowname_col = "set") |>
    gt::tab_header(
      title = "AF Summary Stats"
    ) |>
    gt::tab_stubhead(label = "Set") |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_stub(rows = TRUE)
    ) |>
    gt::tab_options(table.align = "left")

  af_plot <-
    ggplot2::ggplot(data = af_both, ggplot2::aes(.data$af)) +
    ggplot2::geom_histogram(stat = "bin", binwidth = 0.01, fill = "#008080") +
    ggplot2::facet_wrap(~ .data$set, scales = "free_y", drop = FALSE) +
    ggplot2::scale_x_continuous(
      name = "Allele Frequency",
      breaks = seq(0, 1, by = 0.1),
      limits = c(0, 1), expand = c(0, 0)
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, vjust = 1, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    ) +
    ggplot2::labs(title = "AF count distribution")

  list(
    af_stats_gt = af_stats_gt,
    af_plot = af_plot
  )
}
