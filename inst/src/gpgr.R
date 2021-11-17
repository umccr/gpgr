#!/usr/bin/env Rscript --vanilla

suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(glue))
suppressPackageStartupMessages(require(gpgr))


p <- argparse::ArgumentParser(description = "GPG Reporting", prog = "gpgr")
subparser_name <- "subparser_name"
subp <- p$add_subparsers(help = "sub-command help", dest = subparser_name)

#--- HRDetect ---#
hrdetect <- subp$add_parser("hrdetect", help = "hrdetect help")
hrdetect$add_argument("--sample", help = "Sample name.", required = TRUE)
hrdetect$add_argument("--snv", help = "Input SNV (VCF format).", required = TRUE)
hrdetect$add_argument("--sv", help = "Input SV (VCF format).", required = TRUE)
hrdetect$add_argument("--cnv", help = "Input CNV (TSV format).", required = TRUE)
hrdetect$add_argument("--out", help = "Output file ['hrdetect.json.gz'].", default = "hrdetect.json.gz")

#--- CHORD ---#
chord <- subp$add_parser("chord", help = "chord help")
chord$add_argument("--sample", help = "Sample name.", required = TRUE)
chord$add_argument("--snv", help = "Input SNV (VCF format).", required = TRUE)
chord$add_argument("--sv", help = "Input SV (VCF format).", required = TRUE)
chord$add_argument("--out", help = "Output file ['chord.json.gz']", default = "chord.json.gz")

args <- p$parse_args()
if (length(args$subparser_name) == 0) {
  p$print_help()
} else if (args$subparser_name == "hrdetect") {
  # print(c("You've called HRDetect. Here are the arguments: ", args))
  gpgr::hrdetect_run(
    nm = args$sample, snvindel_vcf = args$snv, sv_vcf = args$sv,
    cnv_tsv = args$cnv, outpath = args$out
  )
} else if (args$subparser_name == "chord") {
  # print(c("You've called CHORD. Here are the arguments: ", args))
  gpgr::chord_run(
    vcf.snv = args$snv, vcf.sv = args$sv,
    sample.name = args$sample, outpath = args$out
  )
} else {
  stop("NO IDEA HOW IT GOT TO THIS...")
}

# hrdetect$print_help()
# chord$print_help()
