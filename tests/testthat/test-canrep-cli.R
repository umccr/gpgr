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

test_that("canrep parses required args without error", {
  p <- canrep_parser()
  args <- p$parse_args(required_args())
  expect_equal(args$tumor_name, "x")
})

test_that("canrep does not accept --dragen_hrd", {
  p <- canrep_parser()
  expect_error(
    p$parse_args(c(required_args(), "--dragen_hrd", "sample.hrdscore.csv")),
    regexp = NULL
  )
})
