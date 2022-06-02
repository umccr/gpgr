#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(cli))
suppressPackageStartupMessages(require(gpgr))
suppressPackageStartupMessages(require(glue))

p <- argparse::ArgumentParser(description = "GPG Reporting", prog = "gpgr")
subparser_name <- "subparser_name"
subp <- p$add_subparsers(help = "sub-command help", dest = subparser_name)

source(system.file("cli/linx.R", package = "gpgr"))
source(system.file("cli/canrep.R", package = "gpgr"))

linx_add_args(subp)
canrep_add_args(subp)

args <- p$parse_args()
if (length(args$subparser_name) == 0) {
  p$print_help()
} else if (args$subparser_name == "linx") {
  linx_parse_args(args)
} else if (args$subparser_name == "canrep") {
  canrep_parse_args(args)
} else {
  stop("Need to specify 'linx' or 'canrep' in the cli...")
}
