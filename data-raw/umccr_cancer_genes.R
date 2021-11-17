# UMCCR Cancer Genes
require(here)
require(readr)
require(dplyr)

tmp <- tempfile()
download.file(
  url = "https://raw.githubusercontent.com/umccr/genes/893a655801ce92715f05517b5052e4e81904e870/panels/umccr_2019-03-20.tsv",
  destfile = tmp
)

readr::read_tsv(tmp) |>
  dplyr::select(symbol, tumorsuppressor, oncogene) |>
  readr::write_tsv(here("inst/extdata/ref/umccr_cancer_genes_2019-03-20.tsv"))

tmp2 <- tempfile()
download.file(
  "https://raw.githubusercontent.com/umccr/workflows/8d06a16a0199ccf94b666f9ec027efce8af1110b/genes/cancer_genes/umccr_cancer_genes.hg38.coding.bed",
  destfile = tmp2
)

readr::read_tsv(tmp2, col_names = c("chr", "start", "end", "symbol"), col_types = "ciic") |>
  readr::write_tsv(here("inst/extdata/ref/umccr_cancer_genes_hg38_coding.bed"))
