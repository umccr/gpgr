# gpgr - Genomics Platform Group Reporting

Contains reports and functions used in the Genomics Platform Group
at the University of Melbourne Centre for Cancer Research.

- [gpgr - Genomics Platform Group Reporting](#gpgr---genomics-platform-group-reporting)
  - [Main modules](#main-modules)
    - [PURPLE](#purple)
    - [HRDetect](#hrdetect)
    - [CHORD](#chord)
    - [MutationalPatterns](#mutationalpatterns)
  - [CLI (v0.2.0)](#cli-v020)

[![Conda install](https://anaconda.org/umccr/r-gpgr/badges/installer/conda.svg)](https://anaconda.org/umccr/r-gpgr)

- See <https://umccr.github.io/gpgr/>

## Main modules

### PURPLE

Read and process output files from
[PURPLE](https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator) - see
vignette at <https://umccr.github.io/gpgr/articles/purple.html>.

### HRDetect

Wraps functionality from the [HRDetect](https://github.com/Nik-Zainal-Group/signature.tools.lib)
framework - see
vignette at <https://umccr.github.io/gpgr/articles/hrdetect.html>.

### CHORD

Wraps functionality from [CHORD](https://github.com/UMCUGenetics/CHORD) - see
vignette at <https://umccr.github.io/gpgr/articles/chord.html>.

### MutationalPatterns

Wraps functionality from [MutationalPatterns](https://github.com/UMCUGenetics/MutationalPatterns) - see
vignette at <https://umccr.github.io/gpgr/articles/mutationalpatterns.html>.

## CLI (v0.2.0)

```text
#------- gpgr -------#
$ inst/src/gpgr.R --help
usage: gpgr [-h] {hrdetect,chord,mutpat} ...

GPG Reporting

positional arguments:
  {hrdetect,chord,mutpat}
                        sub-command help
    hrdetect            hrdetect help
    chord               chord help
    mutpat              mutationalpatterns help

#------- hrdetect -------#
$ inst/src/gpgr.R hrdetect --help
usage: gpgr hrdetect [-h] --sample SAMPLE --snv SNV --sv SV --cnv CNV
                     [--out OUT]

optional arguments:
  --sample SAMPLE  Sample name.
  --snv SNV        Input SNV (VCF format).
  --sv SV          Input SV (VCF format).
  --cnv CNV        Input CNV (TSV format).
  --out OUT        Output file ['hrdetect.json.gz'].

#------- chord -------#
$ inst/src/gpgr.R chord --help
usage: gpgr chord [-h] --sample SAMPLE --snv SNV --sv SV [--out OUT]

optional arguments:
  --sample SAMPLE  Sample name.
  --snv SNV        Input SNV (VCF format).
  --sv SV          Input SV (VCF format).
  --out OUT        Output file ['chord.json.gz']

#------- mutpat -------#
$ inst/src/gpgr.R mutpat --help
usage: gpgr mutpat [-h] --sample SAMPLE --snv SNV --outdir OUTDIR

optional arguments:
  --sample SAMPLE  Sample name.
  --snv SNV        Input SNV file (VCF format).
  --outdir OUTDIR  Output directory to write results to
```
