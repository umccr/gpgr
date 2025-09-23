canrep_add_args <- function(subp) {
  canrep <- subp$add_parser("canrep", help = "UMCCR Cancer Report.")
  canrep$add_argument("--af_global", help = "Path to `af_tumor.txt` file.", required = TRUE)
  canrep$add_argument("--af_keygenes", help = "Path to `af_tumor_keygenes.txt` file.", required = TRUE)
  canrep$add_argument("--batch_name", help = "Name of batch sample.", required = TRUE)
  canrep$add_argument("--conda_list", help = "Path to `conda_pkg_list.txt` file.")
  canrep$add_argument("--img_dir", help = "Path to directory containing PURPLE plots.", required = TRUE)
  canrep$add_argument("--key_genes", help = "Path to UMCCR cancer gene file.", required = TRUE)
  canrep$add_argument("--oncokb_genes", help = "Path to OncoKB database file.", required = TRUE)
  canrep$add_argument("--somatic_snv_vcf", help = "Path to `somatic-PASS.vcf.gz` SNV VCF.", required = TRUE)
  canrep$add_argument("--somatic_snv_summary", help = "Path to `somatic_snv_summary.json` JSON.", required = TRUE)
  canrep$add_argument("--somatic_sv_tsv", help = "Path to `manta.tsv` TSV file.", required = TRUE)
  canrep$add_argument("--somatic_sv_vcf", help = "Path to `manta.vcf.gz` VCF file.", required = TRUE)
  canrep$add_argument("--purple_som_gene_cnv", help = "Path to `purple.cnv.gene.tsv`.", required = TRUE)
  canrep$add_argument("--purple_som_cnv_ann", help = "Path to annotated and prioritised `purple.cnv.somatic.tsv`.", required = TRUE)
  canrep$add_argument("--purple_som_cnv", help = "Path to `purple.cnv.somatic.tsv`.", required = TRUE)
  canrep$add_argument("--purple_purity", help = "Path to `purple.purity.tsv`.", required = TRUE)
  canrep$add_argument("--purple_qc", help = "Path to `purple.qc`.", required = TRUE)
  canrep$add_argument("--purple_som_snv_vcf", help = "Path to `purple.somatic.vcf.gz`.", required = TRUE)
  canrep$add_argument("--virusbreakend_tsv", help = "Path to VIRUSBreakend summary file.", required = TRUE)
  canrep$add_argument("--virusbreakend_vcf", help = "Path to VIRUSBreakend VCF file.", required = TRUE)
  canrep$add_argument("--dragen_hrd", help = "Path to DRAGEN HRD file", required = FALSE)
  canrep$add_argument("--bcftools_stats", help = "Path to bcftools stats file", required = TRUE)
  canrep$add_argument("--out_file", help = "Path to output HTML file (needs '.html' suffix) [def: {tumor_name}_cancer_report.html].")
  canrep$add_argument("--quiet", help = "Suppress log printing during rendering.", action = "store_true")
  canrep$add_argument("--result_outdir", help = "Path to directory to write tidy JSON/TSV results.", required = TRUE)
  canrep$add_argument("--tumor_name", help = "Name of tumor sample.", required = TRUE)
}

canrep_parse_args <- function(args) {
  cli::cli_h1("Start rendering UMCCR Cancer Report!")
  if (is.null(args$dragen_hrd)) {
    cli::cli_warn("No DRAGEN HRD file supplied; the report will show DRAGEN HRD scores as missing.")
  }
  res <- gpgr::cancer_rmd(
    af_global = args$af_global,
    af_keygenes = args$af_keygenes,
    batch_name = args$batch_name,
    conda_list = args$conda_list,
    img_dir = args$img_dir,
    key_genes = args$key_genes,
    oncokb_genes = args$oncokb_genes,
    somatic_snv_vcf = args$somatic_snv_vcf,
    somatic_snv_summary = args$somatic_snv_summary,
    somatic_sv_tsv = args$somatic_sv_tsv,
    somatic_sv_vcf = args$somatic_sv_vcf,
    purple_som_gene_cnv = args$purple_som_gene_cnv,
    purple_som_cnv_ann = args$purple_som_cnv_ann,
    purple_som_cnv = args$purple_som_cnv,
    purple_purity = args$purple_purity,
    purple_qc = args$purple_qc,
    purple_som_snv_vcf = args$purple_som_snv_vcf,
    virusbreakend_tsv = args$virusbreakend_tsv,
    virusbreakend_vcf = args$virusbreakend_vcf,
    dragen_hrd = args$dragen_hrd,
    bcftools_stats = args$bcftools_stats,
    out_file = args$out_file,
    quiet = args$quiet,
    result_outdir = args$result_outdir,
    tumor_name = args$tumor_name
  )
  cli::cli_h1("Finished rendering UMCCR Cancer Report!")
  cli::cli_alert_info("Path to HTML output:\n{res}")
}
