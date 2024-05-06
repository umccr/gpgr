#' Generate UMCCR Cancer Report
#'
#' Generates a UMCCR Cancer Report. It does so with the following steps:
#' 1. move the img_dir into 'tmp/img_dir'
#' 2. copy the rmd into 'tmp/cancer_report.Rmd'
#' 3. render the rmd inside 'tmp/'
#' 4. return the path to the output HTML
#'
#' @param af_global Path to `af_tumor.txt` file.
#' @param af_keygenes Path to `af_tumor_keygenes.txt` file.
#' @param batch_name Name of batch sample.
#' @param conda_list Path to `conda_pkg_list.txt` file.
#' @param img_dir Path to directory containing PURPLE plots.
#' @param key_genes Path to UMCCR cancer gene file.
#' @param oncokb_genes Path to OncoKB database file.
#' @param virusbreakend_tsv Path to VIRUSBreakend summary file.
#' @param virusbreakend_vcf Path to VIRUSBreakend VCF file.
#' @param purple_purity Path to `purple.purity.tsv`.
#' @param purple_qc Path to `purple.qc`.
#' @param purple_som_cnv Path to `purple.cnv.somatic.tsv`.
#' @param purple_som_cnv_ann Path to annotated and prioritised `purple.cnv.somatic.tsv`.
#' @param purple_som_gene_cnv Path to `purple.cnv.gene.tsv`.
#' @param purple_som_snv_vcf Path to `purple.somatic.vcf.gz`.
#' @param result_outdir Path to directory to write tidy JSON/TSV results.
#' @param somatic_snv_vcf Path to `somatic-PASS.vcf.gz` SNV VCF.
#' @param somatic_snv_summary Path to `somatic_snv_summary.json` JSON.
#' @param somatic_sv_tsv Path to `manta.tsv` TSV file.
#' @param somatic_sv_vcf Path to `manta.vcf.gz` VCF file.
#' @param tumor_name Name of tumor sample.
#' @param out_file Path to output HTML file (needs '.html' suffix) (def: `{tumor_name}_cancer_report.html`).
#' @param quiet Suppress log printing during rendering.
#' @param bcftools_stats Path to `bcftools_stats.txt` file.
#' @param dragen_hrd Path to DRAGEN `hrdscore.csv` file.
#'
#' @return Path to rendered HTML report.
#' @export
cancer_rmd <- function(af_global, af_keygenes, batch_name, bcftools_stats, conda_list, dragen_hrd, img_dir, key_genes,
                       oncokb_genes, virusbreakend_tsv, virusbreakend_vcf, purple_purity, purple_qc,
                       purple_som_cnv_ann, purple_som_cnv, purple_som_gene_cnv, purple_som_snv_vcf,
                       somatic_snv_vcf, somatic_snv_summary, somatic_sv_tsv, somatic_sv_vcf,
                       result_outdir, tumor_name, out_file = NULL, quiet = FALSE) {
  assertthat::assert_that(
    dir.exists(img_dir),
    quiet %in% c(FALSE, TRUE)
  )
  if (!is.null(out_file)) {
    assertthat::assert_that(
      is.character(out_file),
      tools::file_ext(out_file) == "html"
    )
  } else {
    out_file <- glue::glue("{tumor_name}_cancer_report.html")
  }
  tmp_dir <- tempdir()
  # R's file.copy('foo', 'bar/baz') copies 'foo' to 'bar/baz/foo'
  rmd_dir <- system.file("rmd/umccrise", package = "gpgr")
  img_dir_b <- basename(img_dir)
  rmd_dir_b <- basename(rmd_dir)
  cpdir(rmd_dir, tmp_dir) # /path/to/rmd/umccrise -> /tmp/umccrise
  cpdir(img_dir, file.path(tmp_dir, rmd_dir_b)) # /path/to/um/img -> /tmp/umccrise/img
  rmd_file <- file.path(tmp_dir, "umccrise", "cancer_report.Rmd")
  out_dir <- dirname(out_file)
  mkdir(out_dir)
  pars <- list(
    af_global = af_global,
    af_keygenes = af_keygenes,
    batch_name = batch_name,
    bcftools_stats = bcftools_stats,
    conda_list = conda_list,
    dragen_hrd = dragen_hrd,
    img_dir = img_dir_b,
    key_genes = key_genes,
    oncokb_genes = oncokb_genes,
    somatic_snv_vcf = somatic_snv_vcf,
    somatic_snv_summary = somatic_snv_summary,
    somatic_sv_tsv = somatic_sv_tsv,
    somatic_sv_vcf = somatic_sv_vcf,
    purple_som_gene_cnv = purple_som_gene_cnv,
    purple_som_cnv_ann = purple_som_cnv_ann,
    purple_som_cnv = purple_som_cnv,
    purple_purity = purple_purity,
    purple_qc = purple_qc,
    purple_som_snv_vcf = purple_som_snv_vcf,
    virusbreakend_tsv = virusbreakend_tsv,
    virusbreakend_vcf = virusbreakend_vcf,
    result_outdir = result_outdir,
    tumor_name = tumor_name
  )

  # suppress DT large size warning
  options(DT.warn.size = FALSE)
  rmarkdown::render(
    input = rmd_file,
    params = pars,
    output_dir = out_dir,
    output_file = I(out_file),
    run_pandoc = TRUE,
    quiet = quiet
  )
}
