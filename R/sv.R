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

