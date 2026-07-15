# Tests for pcgr_format_categories() logic in cancer_report.Rmd (sash #52).
# Logic lives inline in the Rmd (can't be sourced directly), so it's
# reimplemented here to get coverage — same pattern as test-hypermutated.R.

pcgr_format_categories <- function(cats_str) {
  cats_str <- paste(cats_str, collapse = ",")
  if (!nzchar(cats_str)) return("")
  tier_labels   <- c("N"="non-coding", "1"="tier 1", "2"="tier 2", "3"="tier 3", "4"="tier 4")
  impact_labels <- c("intergenic"="intergenic", "intronic"="intronic",
                     "downstream"="downstream gene", "upstream"="upstream gene",
                     "impacts_other"="other VEP consequence")
  region_labels <- c("none"="outside GIAB/difficult regions", "difficult"="difficult region",
                     "giab_conf"="GIAB confident region")
  safe_lookup <- function(map, key) {
    val <- unname(map[key])[1]
    if (is.na(val)) key else val
  }
  entries <- trimws(strsplit(cats_str, ",")[[1]])
  entries <- entries[nzchar(entries)]
  lines <- vapply(entries, function(e) {
    parts <- strsplit(e, "\\|")[[1]]
    if (length(parts) != 3) return(NA_character_)
    rc <- strsplit(parts[[3]], ":")[[1]]
    if (length(rc) != 2 || is.na(suppressWarnings(as.integer(rc[[2]])))) return(NA_character_)
    tier   <- parts[[1]]
    impact <- parts[[2]]
    region <- rc[[1]]
    count  <- format(as.integer(rc[[2]]), big.mark = ",", trim = TRUE)
    tl <- safe_lookup(tier_labels, tier)
    il <- safe_lookup(impact_labels, impact)
    rl <- safe_lookup(region_labels, region)
    glue::glue("- {tl} / {il} / {rl}: {count} variants")
  }, character(1), USE.NAMES = FALSE)
  paste(lines[!is.na(lines)], collapse = "\n")
}

test_that("pcgr_format_categories parses a single category", {
  out <- pcgr_format_categories("N|intronic|difficult:80321")
  expect_equal(out, "- non-coding / intronic / difficult region: 80,321 variants")
})

test_that("pcgr_format_categories parses multiple comma-separated categories", {
  out <- pcgr_format_categories("N|intronic|difficult:80321,2|impacts_other|giab_conf:5000")
  expect_equal(
    out,
    paste(
      "- non-coding / intronic / difficult region: 80,321 variants",
      "- tier 2 / other VEP consequence / GIAB confident region: 5,000 variants",
      sep = "\n"
    )
  )
})

test_that("pcgr_format_categories joins multi-element vector input with commas, not fused together", {
  # bolt may emit this as a character vector rather than one pre-joined string
  out <- pcgr_format_categories(c("N|intronic|difficult:80321", "2|impacts_other|giab_conf:5000"))
  expect_equal(
    out,
    paste(
      "- non-coding / intronic / difficult region: 80,321 variants",
      "- tier 2 / other VEP consequence / GIAB confident region: 5,000 variants",
      sep = "\n"
    )
  )
})

test_that("pcgr_format_categories returns empty string for empty input", {
  expect_equal(pcgr_format_categories(""), "")
  expect_equal(pcgr_format_categories(character(0)), "")
})

test_that("pcgr_format_categories skips a malformed entry instead of erroring", {
  # missing the region:count segment on the second entry
  out <- pcgr_format_categories("N|intronic|difficult:80321,malformed_entry")
  expect_equal(out, "- non-coding / intronic / difficult region: 80,321 variants")
})

test_that("pcgr_format_categories skips an entry with a non-numeric count instead of erroring", {
  out <- pcgr_format_categories("N|intronic|difficult:not_a_number")
  expect_equal(out, "")
})

test_that("pcgr_format_categories falls back to the raw code for unknown labels", {
  out <- pcgr_format_categories("9|unknown_impact|unknown_region:100")
  expect_equal(out, "- 9 / unknown_impact / unknown_region: 100 variants")
})
