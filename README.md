
- [📚 gpgr - Genomics Platform Group
  Reporting](#-gpgr---genomics-platform-group-reporting)
  - [Installation](#installation)
  - [Main modules](#main-modules)
    - [🔮 PURPLE](#id_-purple)
    - [🐍 umccrise](#id_-umccrise)
  - [🥳 Developers](#id_-developers)
  - [💻 CLI](#id_-cli)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 📚 gpgr - Genomics Platform Group Reporting

Contains reports and functions used in the Genomics Platform Group at
the University of Melbourne Centre for Cancer Research.

- Docs: <https://umccr.github.io/gpgr/>

[![Conda
install](https://anaconda.org/umccr/r-gpgr/badges/version.svg)](https://anaconda.org/umccr/r-gpgr)
[![Conda
install](https://anaconda.org/umccr/r-gpgr/badges/latest_release_date.svg)](https://anaconda.org/umccr/r-gpgr)

## Installation

``` r
remotes::install_github("umccr/gpgr")
```

- Or if used inside a conda environment:

``` bash
conda install r-gpgr -c umccr -c conda-forge -c bioconda
```

## Main modules

### 🔮 PURPLE

- Read and process output files from the `PURPLE` purity/copy number
  estimator tool from the Hartwig Medical Foundation
  (<https://github.com/hartwigmedical/hmftools/tree/master/purple>). See
  vignette at <https://umccr.github.io/gpgr/articles/purple.html>.

### 🐍 umccrise

- Generate a HTML report with results from the `umccrise` DRAGEN
  tumor/normal post-processing workflow from UMCCR -
  <https://github.com/umccr/umccrise>. See the [CLI](#cli) section below
  for options.

## 🥳 Developers

See <https://umccr.github.io/gpgr/articles/devnotes.html> for developer
notes.

## 💻 CLI

A `gpgr` command line interface is available for convenience.

- If you’re using the conda package, the `gpgr.R` command will already
  be set up inside an activated conda environment.
- If you’re *not* using the conda package, you need to export the
  `gpgr/inst/cli/` directory to your `PATH` in order to use `gpgr.R`.

``` bash
gpgr_cli=$(Rscript -e 'x = system.file("cli", package = "gpgr"); cat(x, "\n")' | xargs)
export PATH="${gpgr_cli}:${PATH}"
```

    gpgr.R --version
    gpgr.R 2.0.0

    #-----------------------------------#
    gpgr.R --help
    usage: gpgr.R [-h] [-v] {canrep} ...

    UMCCR Genomics Platform Group Reporting

    positional arguments:
      {canrep}       sub-command help
        canrep       UMCCR Cancer Report.

    options:
      -h, --help     show this help message and exit
      -v, --version  show program's version number and exit

    #-----------------------------------#
    #------- Cancer Report -------#
    gpgr.R canrep --help
    usage: gpgr.R canrep [-h] --af_global AF_GLOBAL --af_keygenes AF_KEYGENES
                         --batch_name BATCH_NAME [--conda_list CONDA_LIST]
                         --img_dir IMG_DIR --key_genes KEY_GENES --oncokb_genes
                         ONCOKB_GENES --somatic_snv_vcf SOMATIC_SNV_VCF
                         --somatic_snv_summary SOMATIC_SNV_SUMMARY
                         --somatic_sv_tsv SOMATIC_SV_TSV --somatic_sv_vcf
                         SOMATIC_SV_VCF --purple_som_gene_cnv PURPLE_SOM_GENE_CNV
                         --purple_som_cnv_ann PURPLE_SOM_CNV_ANN --purple_som_cnv
                         PURPLE_SOM_CNV --purple_purity PURPLE_PURITY --purple_qc
                         PURPLE_QC --purple_som_snv_vcf PURPLE_SOM_SNV_VCF
                         --virusbreakend_tsv VIRUSBREAKEND_TSV --virusbreakend_vcf
                         VIRUSBREAKEND_VCF --dragen_hrd DRAGEN_HRD
                         --bcftools_stats BCFTOOLS_STATS [--out_file OUT_FILE]
                         [--quiet] --result_outdir RESULT_OUTDIR --tumor_name
                         TUMOR_NAME

    options:
      -h, --help            show this help message and exit
      --af_global AF_GLOBAL
                            Path to `af_tumor.txt` file.
      --af_keygenes AF_KEYGENES
                            Path to `af_tumor_keygenes.txt` file.
      --batch_name BATCH_NAME
                            Name of batch sample.
      --conda_list CONDA_LIST
                            Path to `conda_pkg_list.txt` file.
      --img_dir IMG_DIR     Path to directory containing PURPLE plots.
      --key_genes KEY_GENES
                            Path to UMCCR cancer gene file.
      --oncokb_genes ONCOKB_GENES
                            Path to OncoKB database file.
      --somatic_snv_vcf SOMATIC_SNV_VCF
                            Path to `somatic-PASS.vcf.gz` SNV VCF.
      --somatic_snv_summary SOMATIC_SNV_SUMMARY
                            Path to `somatic_snv_summary.json` JSON.
      --somatic_sv_tsv SOMATIC_SV_TSV
                            Path to `manta.tsv` TSV file.
      --somatic_sv_vcf SOMATIC_SV_VCF
                            Path to `manta.vcf.gz` VCF file.
      --purple_som_gene_cnv PURPLE_SOM_GENE_CNV
                            Path to `purple.cnv.gene.tsv`.
      --purple_som_cnv_ann PURPLE_SOM_CNV_ANN
                            Path to annotated and prioritised
                            `purple.cnv.somatic.tsv`.
      --purple_som_cnv PURPLE_SOM_CNV
                            Path to `purple.cnv.somatic.tsv`.
      --purple_purity PURPLE_PURITY
                            Path to `purple.purity.tsv`.
      --purple_qc PURPLE_QC
                            Path to `purple.qc`.
      --purple_som_snv_vcf PURPLE_SOM_SNV_VCF
                            Path to `purple.somatic.vcf.gz`.
      --virusbreakend_tsv VIRUSBREAKEND_TSV
                            Path to VIRUSBreakend summary file.
      --virusbreakend_vcf VIRUSBREAKEND_VCF
                            Path to VIRUSBreakend VCF file.
      --dragen_hrd DRAGEN_HRD
                            Path to DRAGEN HRD file
      --bcftools_stats BCFTOOLS_STATS
                            Path to bcftools stats file
      --out_file OUT_FILE   Path to output HTML file (needs '.html' suffix) [def:
                            {tumor_name}_cancer_report.html].
      --quiet               Suppress log printing during rendering.
      --result_outdir RESULT_OUTDIR
                            Path to directory to write tidy JSON/TSV results.
      --tumor_name TUMOR_NAME
                            Name of tumor sample.

    #-----------------------------------#
