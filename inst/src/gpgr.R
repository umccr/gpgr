#!/usr/bin/env Rscript --vanilla

suppressPackageStartupMessages(require(argparse))
suppressPackageStartupMessages(require(glue))


p <- argparse::ArgumentParser(description = "GPG Reporting", prog = "gpgr")
subparser_name <- "subparser_name"
subp <- p$add_subparsers(help = "sub-command help", dest=subparser_name)

#--- HRDetect ---#
hrdetect <- subp$add_parser("hrdetect", help="hrdetect help")
hrdetect$add_argument("--sample", help="Sample name.", required=TRUE)
hrdetect$add_argument("--snv", help="Input SNV (VCF format).", required=TRUE)
hrdetect$add_argument("--sv", help="Input SV (VCF format).", required=TRUE)
hrdetect$add_argument("--cnv", help="Input CNV (TSV format).", required=TRUE)
hrdetect$add_argument("--outdir", help="Output directory.", default="hrdetect_output")

#--- CHORD ---#
chord <- subp$add_parser("chord", help="chord help")
chord$add_argument("--sample", help="Sample name.", required=TRUE)
chord$add_argument("--snv", help="Input SNV (VCF format).", required=TRUE)
chord$add_argument("--sv", help="Input SV (VCF format).", required=TRUE)
chord$add_argument("--outdir", help="Output directory.", default="chord_output", required=TRUE)

args <- p$parse_args()
if (length(args$subparser_name) == 0) {
  p$print_help()
} else if (args$subparser_name == "hrdetect") {
  # call hrdetect wrapper
  print(c("You've called HRDetect. Here are the arguments: ", args))
  #gpgr::hrdetect_run(nm=args$sample, snvindel_vcf=args$snv, sv_vcf=args$sv, cnv_tsv=args$cnv, genome="hg38", snvoutdir)

} else if (args$subparser_name == "chord") {
  # call chord wrapper
  print(c("You've called CHORD. Here are the arguments: ", args))
  #gpgr::chord_run(vcf.snv=args$snv, vcf.sv=args$sv, sample.name=args$sample) # a bit slower
} else {
  stop("NO IDEA HOW IT GOT TO THIS...")
}

# hrdetect$print_help()
# chord$print_help()

