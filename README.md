# gpgr - Genomics Platform Group Reporting

Contains reports and functions used in the Genomics Platform Group
at the University of Melbourne Centre for Cancer Research.

[![Conda install](https://anaconda.org/umccr/r-gpgr/badges/installer/conda.svg)](https://anaconda.org/umccr/r-gpgr)

- See <https://umccr.github.io/gpgr/>

## Installation

```r
remotes::install_github("umccr/gpgr")
```

- Or if used inside a conda environment:

```bash
conda install r-gpgr -c umccr -c conda-forge -c bioconda
```

## Main modules

### PURPLE

Read and process output files from
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) - see
vignette at <https://umccr.github.io/gpgr/articles/purple.html>.

### LINX

Read and process output files from
[LINX](https://github.com/hartwigmedical/hmftools/tree/master/linx) - see
vignette at <https://umccr.github.io/gpgr/articles/linx.html>.

## CLI

```bash
# alias to cli
gpgr_cli=$(Rscript -e 'x = system.file("cli", package = "gpgr"); cat(x, "\n")')
alias gpgr=$(which $(echo $gpgr_cli))
```

```text
$ gpgr --help
usage: gpgr [-h] {linx} ...

GPG Reporting

positional arguments:
  {linx}      sub-command help
    linx      LINX HTML report.
```

```text
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
```
