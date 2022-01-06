library(arrow, include.only = "write_parquet")
library(dplyr, include.only = c("bind_rows", "left_join"))
library(readr, include.only = "write_tsv")
library(RJSONIO, include.only = "fromJSON")
source("functions.R")

j <- "multiqc_data.json"
p <- fromJSON(j)
nm <- "umccr_subj_id"
gen <- parse_gen(p) |>
  remove_control_samples() |>
  bind_rows(.id = nm)
raw <- parse_raw(p) |>
  remove_control_samples() |>
  bind_rows(.id = nm)

d <- dplyr::left_join(gen, raw, by = nm)
# data is in tidy format:
# - each variable has its own column
# - each observation (sample) has its own row
# - each value has its own cell

# write to TSV and parquet format
write_parquet(d, "test.parquet")
write_tsv(d, "test.tsv")
