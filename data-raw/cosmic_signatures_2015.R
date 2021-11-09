require(dplyr)
require(readr)
#' Get COSMIC 2015 Signatures
#'
#' Gets COSMIC signature matrix from
#' http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt
#' and returns it as a matrix for signature refitting.
#'
#' @param prob_file Path to signature probabilities file from
#' <http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt>.
#' If not specified, downloads the file from that link.
#'
#' @return A matrix with 96 rows (one for each somatic mutation type within a specific context)
#' and 30 columns (one for each signature).
#'
get_2015_signatures <- function(prob_file = "http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt") {

  # better be explicit - the sig_probs file has 7 extra empty columns
  ctypes <- paste0(c("ccc", paste0(rep("d", 30), collapse = ""), "ccccccc"), collapse = "")
  cnames <- c("SubstType", "Trinucleotide", "SomMutType", paste0("Sig", 1:30), paste0("foo", 1:7))
  cancer_signatures_2015 <-
    readr::read_tsv(prob_file, col_names = cnames, col_types = ctypes, skip = 1) |>
    dplyr::arrange(.data$SubstType)

  # select only 30 sig columns, 96 mut types
  cancer_signatures_2015 |>
    dplyr::select(4:33) |>
    as.matrix()
}

cosmic_signatures_2015 <- get_2015_signatures()

usethis::use_data(cosmic_signatures_2015)
