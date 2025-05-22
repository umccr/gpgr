#' Parse bcftools stats File
#'
#' Parses bcftools stats `bcftools_stats.txt` file.
#'
#' @param x Path to bcftools stats `bcftools_stats.txt` file.
#' @return A ggplot2 object.
#' @export
bcftools_stats_plot <- function(x = NULL) {
  if (is.null(x)) {
    return(NULL)
  }
  cnames <- c("QUAL_dummy", "id", "qual", "snps", "transi", "transv", "indels")
  ln <- readr::read_lines(x)
  d1 <- ln[grepl("QUAL", ln)]
  # line1 ignore, line2 is colnames, just clean those up
  d <- d1[3:length(d1)] |>
    I() |>
    readr::read_tsv(
      col_names = cnames, col_types = readr::cols(.default = "d", "QUAL_dummy" = "c")
    ) |>
    dplyr::select(-"QUAL_dummy")
  if (nrow(d) == 0) {
    return(NULL)
  }
  d <- d |>
    dplyr::select("qual", "snps", "indels") |>
    tidyr::uncount(.data$snps + .data$indels) |>
    dplyr::select("qual")
  med <- stats::median(d$qual, na.rm = TRUE)
  tot <- nrow(d)
  p <- d |>
    ggplot2::ggplot(ggplot2::aes(x = .data$qual)) +
    ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(.data$density)), binwidth = 4, fill = "lightblue") +
    ggplot2::geom_density(alpha = 0.6) +
    ggplot2::geom_vline(xintercept = med, colour = "blue", linetype = "dashed") +
    ggplot2::scale_x_continuous(n.breaks = 10) +
    ggplot2::annotate(
      geom = "label", x = med + 1, y = +Inf, vjust = 2,
      label = paste0("Median: ", med),
    ) +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(glue::glue("Small variant quality score distribution (total variants: {tot})"))
  p
}


#' Parse DRAGEN HRD File
#'
#' Parses DRAGEN `hrdscore.csv` file.
#'
#' @param x Path to DRAGEN `hrdscore.csv` file.
#'
#' @return Tibble with a single row and the following score columns:
#' - HRD (homologous recombination deficiency)
#' - LOH (loss of heterozygosity)
#' - TAI (telomeric allelic imbalance)
#' - LST (large-scale state transitions)
#'
#' If no input is provided, the fields are NA.
#' @export
dragen_hrd <- function(x = NULL) {
  res <- tibble::tibble(Sample = NA, LOH = NA, TAI = NA, LST = NA, HRD = NA)
  if (!is.null(x)) {
    res <- x |>
      readr::read_csv(col_types = readr::cols(.default = "c"))
    names(res) <- sub("_Score", "", names(res))
  }
  res
}

#' HRDetect and CHORD Summary Table
#'
#' Return a summary table with HRDetect and CHORD results.
#'
#' @param chord_res Result from running [sigrap::chord_run()].
#' @param hrdetect_res Result from running [sigrap::hrdetect_run()].
#' @param dragen_res Result from running [gpgr::dragen_hrd()].
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
#' dragen_res <- gpgr::dragen_hrd("path/to/sample.hrdscore.csv")
#' hrd_results_tabs(hrdetect_res = hrdetect_res, chord_res = chord_res, dragen_res = dragen_res)
#' }
#'
#' @export
hrd_results_tabs <- function(hrdetect_res, chord_res, dragen_res) {
  # sn <- chord_res$prediction[, "sample", drop = T]
  # assertthat::are_equal(hrdetect_res[, "sample", drop = T], sn)
    sn <- hrdetect_res[, "sample", drop = TRUE]

  hrdetect_res_tab <-
    tibble::tibble(
      col = colnames(hrdetect_res),
      val = unlist(hrdetect_res[1, ])
    ) |>
    dplyr::filter(.data$col != "sample")

  # chord_res_tab <-
  #   tibble::tibble(
  #     col = colnames(chord_res$prediction),
  #     val = unlist(chord_res$prediction[1, ])
  #   ) |>
  #   dplyr::filter(col != "sample")

  dragen_res_tab <-
    dragen_res |>
    dplyr::select("HRD", "LOH", "TAI", "LST") |>
    unlist() |>
    tibble::enframe(name = "col", value = "val")

  colnames(dragen_res_tab) <- c("DRAGEN", "results_dragen")
  colnames(hrdetect_res_tab) <- c("HRDetect", "results_hrdetect")
  #colnames(chord_res_tab) <- c("CHORD", "results_chord")

  # idea is to make tibbles and cbind them
  tab1 <-
    dplyr::bind_rows(
      hrdetect_res_tab,
      tibble::tribble(
        ~HRDetect, ~results_hrdetect,
        " ", " "
      )
    )
  # tab2a <-
  #   dplyr::bind_rows(
  #     chord_res_tab[1:7, ],
  #     tibble::tribble(
  #       ~CHORD, ~results_chord,
  #       " ", " ",
  #       " ", " "
  #     )
  #   )
  # tab2b <- chord_res_tab[8:16, ] |>
  #   purrr::set_names(c("CHORD2", "results_chord2"))
  tab3 <- dplyr::bind_rows(
    dragen_res_tab,
    tibble::tibble(DRAGEN = rep(" ", 5), results_dragen = rep(" ", 5))
  )

  hrd_results_tab <- dplyr::bind_cols(tab3, tab1)

  hrd_results_gt <-
    hrd_results_tab |>
    gt::gt() |>
    gt::tab_header(
      title = glue::glue("HRD Results for {sn}")
    ) |>
    gt::tab_spanner(
      label = "DRAGEN",
      id = "id_dragen",
      columns = c("DRAGEN", "results_dragen")
    ) |>
    gt::tab_spanner(
      label = "HRDetect",
      id = "id_hrdetect",
      columns = c("HRDetect", "results_hrdetect")
    ) |>
    # gt::tab_spanner(
    #   label = "CHORD",
    #   id = "id_chord",
    #   columns = c("CHORD", "results_chord", "CHORD2", "results_chord2")
    # ) |>
    gt::cols_label(
      DRAGEN = "",
      results_dragen = "",
      HRDetect = "",
      results_hrdetect = "",
      # CHORD = "",
      # results_chord = "",
      # CHORD2 = "",
      # results_chord2 = ""
    ) |>
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        #columns = c("DRAGEN", "HRDetect", "CHORD", "CHORD2")
        columns = c("DRAGEN", "HRDetect")
      )
    ) |>
    gt::cols_align(
      align = "right",
      columns = c("results_hrdetect", "results_dragen")
    ) |>
    gt::tab_style(
      style = gt::cell_borders(
        sides = "left",
        color = "#BBBBBB",
        weight = gt::px(2.5),
        style = "solid"
      ),
      locations = gt::cells_body(
        #columns = c("CHORD", "HRDetect"),
        columns = c("HRDetect"),
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
    dplyr::select("af") |>
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
