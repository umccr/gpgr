#' Read OncoKB
#'
#' @param x Path to something.
#'
#' @return Vector of genes.
#' @export
read_oncokb <- function(x) {
  readr::read_tsv(x, col_types = readr::cols(.default = "c")) |>
    dplyr::filter(
      .data$`OncoKB Annotated` == "Yes"
    ) |>
    dplyr::pull("Hugo Symbol")
}

#' Get OncoKB Genes From Somewhere
#'
#' @param x Path to something.
#' @param oncokb_genes A tibble of something.
#' @return A vector I think.
#'
#' @export
get_oncokb_genes <- function(x, oncokb_genes) {
  delimiters <- " ,&-"
  delimiter_re <- paste0("[", delimiters, "]")

  oncokb_genes |>
    # Create regexes for each match, utilising delimiters for boundaries. Handles most cases where a gene symbol contains the '-' delimiter
    purrr::map(function(n) paste0("(?<=^|", delimiter_re, ")", n, "(?=", delimiter_re, "|$)")) |>
    # Loop with nm iterations through regex and gene symbols
    purrr::map(function(n) stringr::str_detect(x, n) |> tibble::as_tibble_row(.name_repair = "unique_quiet")) |>
    # Combine as tibble to access dplyr::summarise and compile list of detected OncoKB gene symbols for each effect
    dplyr::bind_rows() |>
    dplyr::summarise(dplyr::across(dplyr::everything(), function(v) {
      paste0(sort(oncokb_genes[v]), collapse = ", ")
    })) |>
    unlist()
}
