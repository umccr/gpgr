#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(cli))
suppressPackageStartupMessages(require(gpgr))
suppressPackageStartupMessages(require(glue))

p <- argparse::ArgumentParser(description = "GPG Reporting", prog = "gpgr")
subparser_name <- "subparser_name"
subp <- p$add_subparsers(help = "sub-command help", dest = subparser_name)

#--- LINX report ---#
linx <- subp$add_parser("linx", help = "LINX HTML report.")
linx$add_argument("--sample", help = "Sample name.", required = TRUE)
linx$add_argument("--plot", help = "Path to LINX plot directory.", required = TRUE)
linx$add_argument("--table", help = "Path to LINX table directory.", required = TRUE)
linx$add_argument("--out", help = "HTML output file name [def: linx_{sample}.html].")
linx$add_argument("--quiet", help = "Suppress log printing during rendering.", action = "store_true")


args <- p$parse_args()
if (length(args$subparser_name) == 0) {
  p$print_help()
} else if (args$subparser_name == "linx") {
  # print(c("You've called linx. Here are the arguments: ", args))
  cli::cli_h1("Start rendering LINX R Markdown report!")
  tab <- normalizePath(args$table)
  pl <- normalizePath(args$plot)
  res <- gpgr::linx_rmd(
    sample = args$sample,
    table_dir = tab,
    plot_dir = pl,
    out_file = args$out,
    quiet = args$quiet
  )
  cli::cli_h1("Finished rendering LINX R Markdown report!")
  cli::cli_alert_info("Path to HTML output:\n{res}")
} else {
  stop("Need to specify linx in the cli...")
}

# linx$print_help()
