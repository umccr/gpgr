# UMCCR Cancer Genes
require(here)
require(readr)
require(dplyr)

tmp <- tempfile()
download.file(
  url = "https://raw.githubusercontent.com/umccr/genes/893a655801ce92715f05517b5052e4e81904e870/panels/umccr_2019-03-20.tsv",
  destfile = tmp)

readr::read_tsv(tmp) %>%
  dplyr::select(symbol, tumorsuppressor, oncogene) %>%
  readr::write_tsv(here("inst/extdata/ref/umccr_cancer_genes_2019-03-20.tsv"))
