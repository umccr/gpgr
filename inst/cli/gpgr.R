#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(cli))
suppressPackageStartupMessages(require(gpgr))
suppressPackageStartupMessages(require(glue))

prog_nm <- "gpgr.R"
gpgr_version <- as.character(packageVersion("gpgr"))
p <- argparse::ArgumentParser(
  description = "UMCCR Genomics Platform Group Reporting",
  prog = prog_nm
)
p$add_argument(
  "-v",
  "--version",
  action = "version",
  version = glue::glue("{prog_nm} {gpgr_version}")
)
subparser_name <- "subparser_name"
subp <- p$add_subparsers(help = "sub-command help", dest = subparser_name)

source(system.file("cli/canrep.R", package = "gpgr"))

canrep_add_args(subp)

args <- p$parse_args()
if (length(args$subparser_name) == 0) {
  p$print_help()
} else if (args$subparser_name == "canrep") {
  canrep_parse_args(args)
} else {
  stop("Need to specify 'canrep' in the cli...")
}
