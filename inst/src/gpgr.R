#!/usr/bin/env Rscript --vanilla

# suppressPackageStartupMessages(require(argparser))
# suppressPackageStartupMessages(require(glue))


p <- argparser::arg_parser(description = "GPG Reporting", name = "gpgr", hide.opts = TRUE)
p <- argparser::add_argument(p, arg = "--output", help = "output file", default = "gpgr_report.html")
p <- argparser::add_argument(p, arg = "--input", help = "input dir (specify multiple times)", nargs = Inf)

argv <- argparser::parse_args(p)
input <- argv$input
output <- argv$output

cat(glue::glue("Input: {input}"), sep = "\n")
cat(glue::glue("Output: {output}"), sep = "\n")
