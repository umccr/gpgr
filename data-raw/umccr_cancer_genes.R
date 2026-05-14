# UMCCR Cancer Genes
require(here, include.only = "here")
require(readr, include.only = "read_tsv")
require(dplyr)
require(glue, include.only = "glue")

version <- "24.03.0"
read_umccr_genes_final_panel <- function(version) {
  repo <- "https://raw.githubusercontent.com/umccr/gene_panels"
  tsv_url <- glue::glue(
    "{repo}/v{version}/somatic_panel/3_final_panel/final_panel.tsv"
  )
  d <- readr::read_tsv(
    tsv_url,
    col_types = readr::cols(.default = "c", "tsgene" = "l", "oncogene" = "l")
  ) |>
    dplyr::select(
      symbol = "hgnc_symbol",
      tumorsuppressor = "tsgene",
      "oncogene"
    )
  d
}
read_umccr_genes_final_panel(version) |>
  readr::write_tsv(here(glue(
    "inst/extdata/ref/umccr_cancer_genes_v{version}.tsv"
  )))
