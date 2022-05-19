#!/usr/bin/env Rscript

library(argparser, include.only = c("arg_parser", "add_argument", "parse_args"))
library(arrow, include.only = "write_parquet")
library(readr, include.only = "write_tsv")
source(system.file("scripts/multiqc/functions.R", package = "gpgr"))

p <- arg_parser(
  description = "Export MultiQC json to tidy data.frame/tsv/parquet", hide.opts = TRUE
)
p <- add_argument(p,
  arg = "--json",
  help = "Path to 'multiqc_data.json'."
)
p <- add_argument(p,
  arg = "--outdir",
  help = "Output directory for results.",
  default = "tidymultiqc"
)
p <- add_argument(p,
  arg = "--name",
  help = "Prefix name for output files."
)
args <- parse_args(p)


stopifnot(
  file.exists(args$json), is.character(args$name), nchar(args$name) > 0,
  is.character(args$outdir), nchar(args$outdir) > 0
)
name <- args$name
outdir <- normalizePath(mkdir(args$outdir))
json <- args$json

# main function
d <- mj2df(json)

## write to TSV and Parquet format
tsv_out <- file.path(outdir, paste(name, "tsv", sep = "."))
parquet_out <- file.path(outdir, paste(name, "parquet", sep = "."))
write_tsv(d, tsv_out)
write_parquet(d, parquet_out)
