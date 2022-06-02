
-   <a href="#gpgr---genomics-platform-group-reporting"
    id="toc-gpgr---genomics-platform-group-reporting">gpgr - Genomics
    Platform Group Reporting</a>
    -   <a href="#installation" id="toc-installation">Installation</a>
    -   <a href="#main-modules" id="toc-main-modules">Main modules</a>
        -   <a href="#purple" id="toc-purple">PURPLE</a>
        -   <a href="#linx" id="toc-linx">LINX</a>
    -   <a href="#cli" id="toc-cli">CLI</a>

<!-- README.md is generated from README.Rmd. Please edit that file -->

# gpgr - Genomics Platform Group Reporting

Contains reports and functions used in the Genomics Platform Group at
the University of Melbourne Centre for Cancer Research.

[![Conda
install](https://anaconda.org/umccr/r-gpgr/badges/installer/conda.svg)](https://anaconda.org/umccr/r-gpgr)

-   See <https://umccr.github.io/gpgr/>

## Installation

``` r
remotes::install_github("umccr/gpgr")
```

-   Or if used inside a conda environment:

``` bash
conda install r-gpgr -c umccr -c conda-forge -c bioconda
```

## Main modules

### PURPLE

Read and process output files from
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purple) -
see vignette at <https://umccr.github.io/gpgr/articles/purple.html>.

### LINX

Read and process output files from
[LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx) -
see vignette at <https://umccr.github.io/gpgr/articles/linx.html>.

## CLI

    $ gpgr --help
    usage: gpgr [-h] {linx,canrep} ...

    GPG Reporting

    positional arguments:
      {linx,canrep}  sub-command help
        linx         UMCCR LINX Report.
        canrep       UMCCR Cancer Report.

    optional arguments:
      -h, --help     show this help message and exit



    -------------------------


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



    -------------------------


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
