#' Session Information Kable
#'
#' Session information kables for vignettes.
#'
#' @param pkgs Vector of R packages to display in the vignette.
#'
#' @return A list with two kables containing information about the platform and
#' the specified packages (`pkgs`).
#' @export
session_info_kable <- function(pkgs) {
  si <- devtools::session_info(include_base = TRUE)
  assertthat::assert_that(all(c("platform", "packages") %in% names(si)))
  si_pl <- unclass(si[["platform"]]) |>
    unlist() |>
    tibble::enframe(name = "name", value = "value")
  si_pkg <- unclass(si[["packages"]]) |>
    dplyr::as_tibble() |>
    dplyr::select(
      "package",
      version = "ondiskversion",
      datestamp = "date",
      "source"
    ) |>
    dplyr::filter(.data$package %in% pkgs)

  list(
    si_pl = knitr::kable(si_pl, caption = "Platform information."),
    si_pkg = knitr::kable(
      si_pkg,
      caption = "Main packages used in this vignette."
    )
  )
}
