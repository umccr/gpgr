#!/usr/bin/env Rscript

require(here)
require(gpgr)
require(glue)
require(purrr)

batch_name <- "SBJ02269__PRJ221202"
tumor_name <- "PRJ221202"
# umccrised_dir <- "/g/data/gx8/projects/diakumis/umccrise/data/SBJ02269"
# key_genes <- "/g/data/gx8/projects/diakumis/conda/envs/umccrise/lib/python3.7/site-packages/ngs_utils/reference_data/key_genes/umccr_cancer_genes.latest.tsv"
umccrised_dir <- here::here("nogit/umccrised_data/SBJ02269")
key_genes <- here::here("nogit/umccrised_data/umccr_cancer_genes.latest.tsv")

params <- list(
  af_global = glue::glue("{umccrised_dir}/work/{batch_name}/cancer_report/afs/af_tumor.txt"),
  af_keygenes = glue::glue("{umccrised_dir}/work/{batch_name}/cancer_report/afs/af_tumor_keygenes.txt"),
  batch_name = batch_name,
  conda_list = glue::glue("{umccrised_dir}/work/{batch_name}/conda_pkg_list.txt"),
  img_dir = glue::glue("{umccrised_dir}/img"),
  key_genes = key_genes,
  oncoviral_breakpoints_tsv = glue::glue("{umccrised_dir}/work/{batch_name}/oncoviruses/oncoviral_breakpoints.tsv"),
  oncoviral_present_viruses = glue::glue("{umccrised_dir}/work/{batch_name}/oncoviruses/present_viruses.txt"),
  out_file = here::here("nogit/umccrised_data/SBJ02269/foo_cancer_report.html"),
  purple_purity = glue::glue("{umccrised_dir}/work/{batch_name}/purple/{batch_name}.purple.purity.tsv"),
  purple_qc = glue::glue("{umccrised_dir}/work/{batch_name}/purple/{batch_name}.purple.qc"),
  purple_som_cnv = glue::glue("{umccrised_dir}/work/{batch_name}/purple/{batch_name}.purple.cnv.somatic.tsv"),
  purple_som_gene_cnv = glue::glue("{umccrised_dir}/work/{batch_name}/purple/{batch_name}.purple.cnv.gene.tsv"),
  purple_som_snv_vcf = glue::glue("{umccrised_dir}/work/{batch_name}/purple/{batch_name}.purple.somatic.vcf.gz"),
  result_outdir = glue::glue("{umccrised_dir}/{batch_name}/cancer_report_tables2"),
  somatic_snv_summary = glue::glue("{umccrised_dir}/work/{batch_name}/cancer_report/somatic_snv_summary.json"),
  somatic_snv_vcf = glue::glue("{umccrised_dir}/{batch_name}/small_variants/{batch_name}-somatic-PASS.vcf.gz"),
  somatic_sv_tsv = glue::glue("{umccrised_dir}/{batch_name}/structural/{batch_name}-manta.tsv"),
  somatic_sv_vcf = glue::glue("{umccrised_dir}/{batch_name}/structural/{batch_name}-manta.vcf.gz"),
  tumor_name = tumor_name
)

# awesomeness
purrr::lift_dl(gpgr::cancer_rmd)(params)
