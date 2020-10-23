#' HRDetect and CHORD Summary Table
#'
#' Return a summary table with HRDetect and CHORD results.
#'
#' @param chord_res Result from running [gpgr::chord_run()].
#' @param hrdetect_res Result from running [gpgr::hrdetect_run()].
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
#' hrdetect_res <- hrdetect_run(nm, snv, sv, cnv, genome, snvoutdir)
#' chord_res <- chord_run(vcf.snv = snv, df.sv = gpgr:::chord_mantavcf2df(sv), sample.name = nm)
#' hrd_results_tabs(hrdetect_res = hrdetect_res, chord_res = chord_res)
#' }
#'
#' @export
hrd_results_tabs <- function(hrdetect_res, chord_res) {

  sn <- chord_res$prediction[, "sample", drop = T]
  assertthat::are_equal(hrdetect_res[, "sample", drop = T], sn)

  hrdetect_res_tab <-
    tibble::tibble(col = colnames(hrdetect_res),
                   val = unlist(hrdetect_res[1, ])) %>%
    dplyr::filter(.data$col != "sample")

  chord_res_tab <-
    tibble::tibble(col = colnames(chord_res$prediction),
                   val = unlist(chord_res$prediction[1, ])) %>%
    dplyr::filter(col != "sample")


  colnames(hrdetect_res_tab) <- c("HRDetect", "results_hrdetect")
  colnames(chord_res_tab) <- c("CHORD", "results_chord")

  tab1 <-
    dplyr::bind_rows(
      hrdetect_res_tab,
      tibble::tribble(
        ~HRDetect, ~results_hrdetect,
        " ", " "))
  tab2 <-
    dplyr::bind_rows(
      chord_res_tab[1:7, ],
      tibble::tribble(
        ~CHORD, ~results_chord,
        " ", " ",
        " ", " "))
  tab3 <- chord_res_tab[8:16, ] %>%
    purrr::set_names(c("CHORD2", "results_chord2"))

  hrd_results_tab <- dplyr::bind_cols(tab1, tab2, tab3)

  hrd_results_gt <-
    hrd_results_tab %>%
    gt::gt() %>%
    gt::tab_header(
      title = glue::glue("HRD Results for {sn}")
    ) %>%
    gt::tab_spanner(
      label = "HRDetect",
      columns =  c("HRDetect", "results_hrdetect")
    ) %>%
    gt::tab_spanner(
      label = "CHORD",
      columns = c("CHORD", "results_chord", "CHORD2", "results_chord2")
    ) %>%
    gt::cols_label(
      HRDetect = "",
      results_hrdetect = "",
      CHORD = "",
      results_chord = "",
      CHORD2 = "",
      results_chord2 = ""
    ) %>%
    gt::tab_style(
      style = list(
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(
        columns = c("HRDetect", "CHORD", "CHORD2")
      )
    ) %>%
    gt::cols_align(
      align = "right",
      columns = c("results_hrdetect")
    ) %>%
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
    )

  list(sample = sn,
       hrd_results_tab = hrd_results_tab,
       hrd_results_gt = hrd_results_gt
  )

}
