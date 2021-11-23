# NEWS

## 0.2.0 (2021-11-23)

- new release after a year in the freezer
- :star: add CLI for CHORD, HRDetect and MutationalPatterns.
  ([pr26](https://github.com/umccr/gpgr/pull/26),
  [pr28](https://github.com/umccr/gpgr/pull/28),
  [issue17](https://github.com/umccr/gpgr/issues/17))
- :star: add 6 new SBS signatures from the March 2021 COSMIC release (v3.2).
  ([issue22](https://github.com/umccr/gpgr/issues/22),
  [pr23](https://github.com/umccr/gpgr/pull/23))
- :star: update to use the base pipe operator (`|>`) introduced in R 4.1.
  ([pr24](https://github.com/umccr/gpgr/pull/24))
- :star: add handy `get_genome_obj` function that returns a BSgenome object, provided
  a 'hgXX' string (useful for MutationalPatterns & CHORD).
  ([pr20](https://github.com/umccr/gpgr/pull/20))
- :star: add handy `pkg_exists` function that checks if the specified R package
  exists on the local system. ([pr20](https://github.com/umccr/gpgr/pull/20))
- :question: add `readr::local_edition(1)` before `readr::read_tsv` due to
  `readr` using `vroom` under the hood (I'm not too comfortable with this,
  but anyway). See <https://www.tidyverse.org/blog/2021/07/readr-2-0-0/#readr-2nd-edition>.
  I think it happens because I'm using the `comment` parameter to read in the
  gzipped VCF. Oh well. ([pr21](https://github.com/umccr/gpgr/pull/21))
- :wrench: use [bumpversion](https://github.com/c4urself/bump2version) for
  updating pkg versions.
  ([issue18](https://github.com/umccr/gpgr/issues/18),
  [pr19](https://github.com/umccr/gpgr/pull/19))
- :wrench: add pre-commit hooks and styler.
  ([issue25](https://github.com/umccr/gpgr/issues/25),
  [pr27](https://github.com/umccr/gpgr/pull/27))

## v0.1.0 (2020-11-19)

- Handle umccrise cancer report functionality:
  - PURPLE CNVs
  - SVs (tables + BND plots)
  - AF (summary + plot)
  - HRD
  - write final tables to tsv + json
  - Support for MutationalPatterns v3.0 (https://github.com/UMCUGenetics/MutationalPatterns)
  - Support for PURPLE v2.51+

## v0.0.6 (2020-10-09)

- Initial release of gpgr. Features include:
  - Support for CHORD (https://github.com/UMCUGenetics/CHORD)
  - Support for HRDetect (https://github.com/Nik-Zainal-Group/signature.tools.lib)
  - Support for PURPLE (https://github.com/hartwigmedical/hmftools/tree/master/purple)
  - R package webpage deployed via GitHub Actions at https://umccr.github.io/gpgr/