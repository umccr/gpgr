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

