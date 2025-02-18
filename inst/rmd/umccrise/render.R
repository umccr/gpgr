#!/usr/bin/env Rscript

require(here, include.only = "here")
require(gpgr)
require(glue, include.only = "glue")
require(purrr)

batch_name <- "SBJ00480_PTC_HCC1395t100pc"
tumor_name <- "PTC_HCC1395t100pc"
umccrised_dir <- here::here("nogit/umccrised_data/seqc_inputs.20231115")
params <- list(
  af_global = glue("{umccrised_dir}/sample_data/af_tumor.txt"),
  af_keygenes = glue("{umccrised_dir}/sample_data/af_tumor_keygenes.txt"),
  batch_name = batch_name,
  bcftools_stats = glue("{umccrised_dir}/sample_data/{tumor_name}.somatic.bcftools_stats.txt"),
  conda_list = NULL,
  dragen_hrd = glue("{umccrised_dir}/sample_data/{tumor_name}.hrdscore.csv"),
  img_dir = glue("{umccrised_dir}/output/img"),
  key_genes = normalizePath("./inst/extdata/ref/umccr_cancer_genes_v24.03.0.tsv"),
  oncokb_genes = glue("{umccrised_dir}/reference_data/oncokb_genes.20231113.tsv"),
  virusbreakend_tsv = glue("{umccrised_dir}/sample_data/virusbreakend/{batch_name}.virusbreakend.vcf.summary.tsv"),
  virusbreakend_vcf = glue("{umccrised_dir}/sample_data/virusbreakend/{batch_name}.virusbreakend.vcf"),
  purple_purity = glue("{umccrised_dir}/sample_data/purple/{tumor_name}.purple.purity.tsv"),
  purple_qc = glue("{umccrised_dir}/sample_data/purple/{tumor_name}.purple.qc"),
  purple_som_cnv_ann = glue("{umccrised_dir}/sample_data/{tumor_name}.cnv.prioritised.tsv"),
  purple_som_cnv = glue("{umccrised_dir}/sample_data/purple/{tumor_name}.purple.cnv.somatic.tsv"),
  purple_som_gene_cnv = glue("{umccrised_dir}/sample_data/purple/{tumor_name}.purple.cnv.gene.tsv"),
  purple_som_snv_vcf = glue("{umccrised_dir}/sample_data/purple/{tumor_name}.purple.somatic.vcf.gz"),
  result_outdir = glue("{umccrised_dir}/output/cancer_report_tables"),
  somatic_snv_vcf = glue("{umccrised_dir}/sample_data/{tumor_name}.pass.vcf.gz"),
  somatic_snv_summary = glue("{umccrised_dir}/sample_data/{tumor_name}.somatic.variant_counts_process.json"),
  somatic_sv_tsv = here::here("inst/extdata/sash/sv.prioritised.tsv"),
  somatic_sv_vcf = here::here("inst/extdata/sash/sv.prioritised.vcf.gz"),
  tumor_name = tumor_name,
  mutpat_dir = NULL,
  hrdetect_file = NULL,
  chord_file = NULL
)

# awesomeness
purrr::lift_dl(gpgr::cancer_rmd)(params)
