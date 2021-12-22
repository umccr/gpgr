#!/usr/bin/env Rscript --vanilla

suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(gpgr))


p <- argparse::ArgumentParser(description = "GPG Reporting", prog = "gpgr")
subparser_name <- "subparser_name"
subp <- p$add_subparsers(help = "sub-command help", dest = subparser_name)

#--- Example ---#
example <- subp$add_parser("example", help = "example help")
example$add_argument("--sample", help = "Sample name.", required = TRUE)
example$add_argument("--snv", help = "Input SNV (VCF format).", required = TRUE)
example$add_argument("--sv", help = "Input SV (VCF format).", required = TRUE)
example$add_argument("--cnv", help = "Input CNV (TSV format).", required = TRUE)
example$add_argument("--out", help = "Output file ['example.txt'].", default = "example.txt")

args <- p$parse_args()
if (length(args$subparser_name) == 0) {
  p$print_help()
} else if (args$subparser_name == "example") {
  # print(c("You've called Example Here are the arguments: ", args))
  gpgr::example_run(
    nm = args$sample, snvindel_vcf = args$snv, sv_vcf = args$sv,
    cnv_tsv = args$cnv, outpath = args$out
  )
} else {
  stop("NO IDEA HOW IT GOT TO THIS...")
}

# example$print_help()
