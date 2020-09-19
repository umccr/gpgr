#' Set Up gpgr Variables
#'
#' Sets up gpgr variables after listing the contents of a umccrise `rmd` directory.
#'
#' @param x Path to umccrise `rmd` directory.
#'
#' @return A list with several elements pointing to umccrise output file names.
#'
#' @examples
#'
#' @testexamples
#'
#' @export
setup_vars <- function(x) {

  assertthat::assert_that(length(x) == 1, is.character(x))


}
