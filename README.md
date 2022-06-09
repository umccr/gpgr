
-   <a href="#-gpgr---genomics-platform-group-reporting"
    id="toc--gpgr---genomics-platform-group-reporting">📚 gpgr - Genomics
    Platform Group Reporting</a>
    -   <a href="#installation" id="toc-installation">Installation</a>
    -   <a href="#main-modules" id="toc-main-modules">Main modules</a>
        -   <a href="#-linx" id="toc--linx">🕸 LINX</a>
        -   <a href="#-purple" id="toc--purple">🔮 PURPLE</a>
        -   <a href="#-umccrise" id="toc--umccrise">🐍 umccrise</a>
    -   <a href="#-cli" id="toc--cli">💻 CLI</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

# 📚 gpgr - Genomics Platform Group Reporting

Contains reports and functions used in the Genomics Platform Group at
the University of Melbourne Centre for Cancer Research.

-   Docs: <https://umccr.github.io/gpgr/>

[![Conda
install](https://anaconda.org/umccr/r-gpgr/badges/installer/conda.svg)](https://anaconda.org/umccr/r-gpgr)

## Installation

``` r
remotes::install_github("umccr/gpgr")
```

-   Or if used inside a conda environment:

``` bash
conda install r-gpgr -c umccr -c conda-forge -c bioconda
```

## Main modules

### 🕸 LINX

-   Generate a HTML report with results from the `LINX` structural
    variant visualisation tool from the Hartwig Medical Foundation
    (<https://github.com/hartwigmedical/hmftools/tree/master/linx>). See
    the [CLI](#cli) section below for options.
-   For useful functions for reading/processing `LINX` results, see the
    vignette at <https://umccr.github.io/gpgr/articles/linx.html>.

### 🔮 PURPLE

-   Read and process output files from the `PURPLE` purity/copy number
    estimator tool from the Hartwig Medical Foundation
    (<https://github.com/hartwigmedical/hmftools/tree/master/purple>).
    See vignette at <https://umccr.github.io/gpgr/articles/purple.html>.

### 🐍 umccrise

-   Generate a HTML report with results from the `umccrise` DRAGEN
    tumor/normal post-processing workflow from UMCCR -
    <https://github.com/umccr/umccrise>. See the [CLI](#cli) section
    below for options.

## 💻 CLI

A `gpgr` command line interface is available for convenience.

-   If you’re using the conda package, you can use the `gpgr.R` command
    inside an activated conda environment.
-   If you want to set it up on your local machine (non-conda), you can
    use the `alias` command chunk below (or just make that script
    available in your `PATH`).

<!-- -->

    $ gpgr --help
    usage: gpgr [-h] [-v] {linx,canrep} ...

    UMCCR Genomics Platform Group Reporting

    positional arguments:
      {linx,canrep}  sub-command help
        linx         UMCCR LINX Report.
        canrep       UMCCR Cancer Report.

    optional arguments:
      -h, --help     show this help message and exit
      -v, --version  show program's version number and exit



    $ gpgr --version
    gpgr 1.2.6



    #------- LINX Report -------#


    $ gpgr linx --help
    usage: gpgr linx [-h] --sample SAMPLE --plot PLOT --table TABLE [--out OUT]
                     [--quiet]

    optional arguments:
      -h, --help       show this help message and exit
      --sample SAMPLE  Sample name.
      --plot PLOT      Path to LINX plot directory.
      --table TABLE    Path to LINX table directory.
      --out OUT        HTML output file name [def: linx_{sample}.html].
      --quiet          Suppress log printing during rendering.



    #------- Cancer Report -------#


    $ gpgr canrep --help
    usage: gpgr canrep [-h] --af_global AF_GLOBAL --af_keygenes AF_KEYGENES
                       --batch_name BATCH_NAME --conda_list CONDA_LIST --img_dir
                       IMG_DIR --key_genes KEY_GENES --somatic_snv_vcf
                       SOMATIC_SNV_VCF --somatic_sv_tsv SOMATIC_SV_TSV
                       --somatic_sv_vcf SOMATIC_SV_VCF --purple_som_gene_cnv
                       PURPLE_SOM_GENE_CNV --purple_som_cnv PURPLE_SOM_CNV
                       --purple_germ_cnv PURPLE_GERM_CNV --purple_purity
                       PURPLE_PURITY --purple_qc PURPLE_QC --purple_som_snv_vcf
                       PURPLE_SOM_SNV_VCF --oncoviral_present_viruses
                       ONCOVIRAL_PRESENT_VIRUSES --oncoviral_breakpoints_tsv
                       ONCOVIRAL_BREAKPOINTS_TSV [--out_file OUT_FILE] [--quiet]
                       --result_outdir RESULT_OUTDIR --tumor_name TUMOR_NAME

    optional arguments:
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
      --somatic_snv_vcf SOMATIC_SNV_VCF
                            Path to `somatic-PASS.vcf.gz` SNV VCF.
      --somatic_sv_tsv SOMATIC_SV_TSV
                            Path to `manta.tsv` TSV file.
      --somatic_sv_vcf SOMATIC_SV_VCF
                            Path to `manta.vcf.gz` VCF file.
      --purple_som_gene_cnv PURPLE_SOM_GENE_CNV
                            Path to `purple.cnv.gene.tsv`.
      --purple_som_cnv PURPLE_SOM_CNV
                            Path to `purple.cnv.somatic.tsv`.
      --purple_germ_cnv PURPLE_GERM_CNV
                            Path to `purple.cnv.germline.tsv`.
      --purple_purity PURPLE_PURITY
                            Path to `purple.purity.tsv`.
      --purple_qc PURPLE_QC
                            Path to `purple.qc`.
      --purple_som_snv_vcf PURPLE_SOM_SNV_VCF
                            Path to `purple.somatic.vcf.gz`.
      --oncoviral_present_viruses ONCOVIRAL_PRESENT_VIRUSES
                            Path to `oncoviruses/present_viruses.txt`.
      --oncoviral_breakpoints_tsv ONCOVIRAL_BREAKPOINTS_TSV
                            Path to `oncoviruses/oncoviral_breakpoints.tsv`.
      --out_file OUT_FILE   Path to output HTML file (needs '.html' suffix) [def:
                            {tumor_name}_cancer_report.html].
      --quiet               Suppress log printing during rendering.
      --result_outdir RESULT_OUTDIR
                            Path to directory to write tidy JSON/TSV results.
      --tumor_name TUMOR_NAME
                            Name of tumor sample.
