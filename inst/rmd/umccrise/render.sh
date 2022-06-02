#!/usr/bin/env bash

set -euo pipefail

gpgr_dir="/Users/pdiakumis/projects/gpgr"
umccrised_dir="${gpgr_dir}/nogit/umccrised_data/SBJ02269"
key_genes="${gpgr_dir}/nogit/umccrised_data/umccr_cancer_genes.latest.tsv"
batch_name="SBJ02269__PRJ221202"
tumor_name="PRJ221202"

# link to cli
gpgr_cli=$(Rscript -e 'x = system.file("cli/gpgr.R", package = "gpgr"); cat(x, "\n")')

$gpgr_cli canrep \
  --af_global "${umccrised_dir}/work/${batch_name}/cancer_report/afs/af_tumor.txt" \
  --af_keygenes "${umccrised_dir}/work/${batch_name}/cancer_report/afs/af_tumor_keygenes.txt" \
  --batch_name "${batch_name}" \
  --conda_list  "${umccrised_dir}/work/${batch_name}/conda_pkg_list.txt" \
  --img_dir "${umccrised_dir}/img" \
  --key_genes  "${key_genes}" \
  --somatic_snv_vcf  "${umccrised_dir}/${batch_name}/small_variants/${batch_name}-somatic-PASS.vcf.gz" \
  --somatic_sv_tsv  "${umccrised_dir}/${batch_name}/structural/${batch_name}-manta.tsv" \
  --somatic_sv_vcf  "${umccrised_dir}/${batch_name}/structural/${batch_name}-manta.vcf.gz" \
  --purple_som_gene_cnv  "${umccrised_dir}/work/${batch_name}/purple/${batch_name}.purple.cnv.gene.tsv" \
  --purple_som_cnv  "${umccrised_dir}/work/${batch_name}/purple/${batch_name}.purple.cnv.somatic.tsv" \
  --purple_germ_cnv  "${umccrised_dir}/work/${batch_name}/purple/${batch_name}.purple.cnv.germline.tsv" \
  --purple_purity  "${umccrised_dir}/work/${batch_name}/purple/${batch_name}.purple.purity.tsv" \
  --purple_qc  "${umccrised_dir}/work/${batch_name}/purple/${batch_name}.purple.qc" \
  --purple_som_snv_vcf  "${umccrised_dir}/work/${batch_name}/purple/${batch_name}.purple.somatic.vcf.gz" \
  --oncoviral_present_viruses  "${umccrised_dir}/work/${batch_name}/oncoviruses/present_viruses.txt" \
  --oncoviral_breakpoints_tsv  "${umccrised_dir}/work/${batch_name}/oncoviruses/oncoviral_breakpoints.tsv" \
  --out_file  "${gpgr_dir}/nogit/umccrised_data/cancer_report_2269.html" \
  --result_outdir  "${umccrised_dir}/${batch_name}/cancer_report_tables2" \
  --tumor_name  tumor_name
