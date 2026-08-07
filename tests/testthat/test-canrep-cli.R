# Tests for canrep CLI argument parsing

canrep_parser <- function() {
  cli_path <- system.file("cli/canrep.R", package = "gpgr")
  source(cli_path, local = TRUE)
  p <- argparse::ArgumentParser()
  sp <- p$add_subparsers(dest = "command")
  canrep_add_args(sp)
  p
}

required_args <- function() {
  c(
    "canrep",
    "--af_global", "x",
    "--af_keygenes", "x",
    "--batch_name", "x",
    "--img_dir", "x",
    "--key_genes", "x",
    "--oncokb_genes", "x",
    "--somatic_snv_vcf", "x",
    "--somatic_snv_summary", "x",
    "--somatic_sv_tsv", "x",
    "--somatic_sv_vcf", "x",
    "--purple_som_gene_cnv", "x",
    "--purple_som_cnv_ann", "x",
    "--purple_som_cnv", "x",
    "--purple_purity", "x",
    "--purple_qc", "x",
    "--purple_som_snv_vcf", "x",
    "--virusbreakend_tsv", "x",
    "--virusbreakend_vcf", "x",
    "--bcftools_stats", "x",
    "--result_outdir", "x",
    "--tumor_name", "x"
  )
}

# Guards against regression where --dragen_hrd was accidentally marked required=TRUE
# (introduced in 563f946, fixed in PR #94 and again in 2.3.1)

test_that("canrep parses without --dragen_hrd (optional)", {
  p <- canrep_parser()
  args <- p$parse_args(required_args())
  expect_null(args$dragen_hrd)
})

test_that("canrep parses with --dragen_hrd when provided", {
  p <- canrep_parser()
  args <- p$parse_args(c(required_args(), "--dragen_hrd", "sample.hrdscore.csv"))
  expect_equal(args$dragen_hrd, "sample.hrdscore.csv")
})
