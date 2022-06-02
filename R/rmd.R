#' Generate UMCCR LINX Report
#'
#' Generates a UMCCR LINX HTML report. It does so with the following steps:
#' 1. recursively copy the table and plot directories into tmp/plot and tmp/table
#' 2. copy the rmd into tmp/linx.Rmd
#' 3. render the rmd inside tmp/
#' 4. return the path to the output HTML
#'
#' @param sample Name of sample.
#' @param table_dir Path to LINX table directory.
#' @param plot_dir Path to LINX plot directory.
#' @param out_file Path to output HTML file (needs '.html' suffix).
#' @param quiet Suppress printing during rendering.
#'
#' @return Path to rendered HTML report.
#' @export
linx_rmd <- function(sample, table_dir, plot_dir, out_file = NULL, quiet = FALSE) {
  assertthat::assert_that(
    is.character(sample), nchar(sample) > 0,
    dir.exists(table_dir), dir.exists(plot_dir),
    quiet %in% c(FALSE, TRUE)
  )
  if (!is.null(out_file)) {
    assertthat::assert_that(
      is.character(out_file),
      tools::file_ext(out_file) == "html"
    )
  } else {
    out_file <- glue::glue("linx_{sample}.html")
  }
  tmp_dir <- tempdir()
  rmd_dir <- system.file("rmd/linx", package = "gpgr")
  cpdir(table_dir, tmp_dir)
  cpdir(plot_dir, tmp_dir)
  cpdir(rmd_dir, tmp_dir)
  rmd_file <- file.path(tmp_dir, "linx", "linx.Rmd")
  out_dir <- dirname(out_file)
  mkdir(out_dir)
  pars <- list(
    table_dir = table_dir,
    plot_dir = plot_dir,
    sample = sample
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

#' Generate UMCCR Cancer Report
#'
#' Generates a UMCCR Cancer Report. It does so with the following steps:
#' 1. recursively copy the img_dir into tmp/img
#' 2. copy the rmd into tmp/cancer_report.Rmd
#' 3. render the rmd inside tmp/
#' 4. return the path to the output HTML
#'
#' @param af_global Path to `af_tumor.txt` file.
#' @param af_keygenes Path to `af_tumor_keygenes.txt` file.
#' @param batch_name Name of batch sample.
#' @param conda_list Path to `conda_pkg_list.txt` file.
#' @param img_dir Path to directory containing PURPLE plots.
#' @param key_genes Path to UMCCR cancer gene file.
#' @param somatic_snv_vcf Path to `somatic-PASS.vcf.gz` SNV VCF.
#' @param somatic_sv_tsv Path to `manta.tsv` TSV file.
#' @param somatic_sv_vcf Path to `manta.vcf.gz` VCF file.
#' @param purple_som_gene_cnv Path to `purple.cnv.gene.tsv`.
#' @param purple_som_cnv Path to `purple.cnv.somatic.tsv`.
#' @param purple_germ_cnv Path to `purple.cnv.germline.tsv`.
#' @param purple_purity Path to `purple.purity.tsv`.
#' @param purple_qc Path to `purple.qc`.
#' @param purple_som_snv_vcf Path to `purple.somatic.vcf.gz`.
#' @param oncoviral_present_viruses Path to `oncoviruses/present_viruses.txt`.
#' @param oncoviral_breakpoints_tsv Path to `oncoviruses/oncoviral_breakpoints.tsv`.
#' @param out_file Path to output HTML file (needs '.html' suffix) (def: `{tumor_name}_cancer_report.html`).
#' @param quiet Suppress log printing during rendering.
#' @param result_outdir Path to directory to write tidy JSON/TSV results.
#' @param tumor_name Name of tumor sample.
#'
#' @return Path to rendered HTML report.
#' @export
cancer_rmd <- function(af_global, af_keygenes, batch_name, conda_list, img_dir, key_genes,
                       somatic_snv_vcf, somatic_sv_tsv, somatic_sv_vcf, purple_som_gene_cnv,
                       purple_som_cnv, purple_germ_cnv, purple_purity, purple_qc,
                       purple_som_snv_vcf, oncoviral_present_viruses, oncoviral_breakpoints_tsv,
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
  rmd_dir <- system.file("rmd/umccrise", package = "gpgr")
  cpdir(img_dir, tmp_dir)
  cpdir(rmd_dir, tmp_dir)
  rmd_file <- file.path(tmp_dir, "umccrise", "cancer_report.Rmd")
  out_dir <- dirname(out_file)
  mkdir(out_dir)
  pars <- list(
    af_global = af_global,
    af_keygenes = af_keygenes,
    batch_name = batch_name,
    conda_list = conda_list,
    img_dir = img_dir,
    key_genes = key_genes,
    somatic_snv_vcf = somatic_snv_vcf,
    somatic_sv_tsv = somatic_sv_tsv,
    somatic_sv_vcf = somatic_sv_vcf,
    purple_som_gene_cnv = purple_som_gene_cnv,
    purple_som_cnv = purple_som_cnv,
    purple_germ_cnv = purple_germ_cnv,
    purple_purity = purple_purity,
    purple_qc = purple_qc,
    purple_som_snv_vcf = purple_som_snv_vcf,
    oncoviral_present_viruses = oncoviral_present_viruses,
    oncoviral_breakpoints_tsv = oncoviral_breakpoints_tsv,
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
