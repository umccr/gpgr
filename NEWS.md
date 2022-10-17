# 1.5.0rc

- add bpi status in canrep summary table

# gpgr 1.4.0

- :bug: Update LINX cluster plot path pattern ([issue50](https://github.com/umccr/gpgr/issues/50),
  [pr51](https://github.com/umccr/gpgr/pull/51)).
- :star: umccrise canrep: add SNV summary ([pr53](https://github.com/umccr/gpgr/pull/53)).
  - add hypermutated status in canrep summary table
- :bug: Handle empty BPI Start/End values ([pr55](https://github.com/umccr/gpgr/pull/55),
  [issue54](https://github.com/umccr/gpgr/issues/54), [umccrise issue88](https://github.com/umccr/umccrise/issues/88)).
- :bug: `{gt}` v0.7.0 requires a unique `id` for each `tab_spanner` ([issue57](https://github.com/umccr/gpgr/issues/57),
  [pr58](https://github.com/umccr/gpgr/pull/58)).

# gpgr 1.3.0 (2022-06-10)

- :star: Add Docker support ([pr49](https://github.com/umccr/gpgr/pull/49)).
- :star: Add conda-lock support ([pr49](https://github.com/umccr/gpgr/pull/49)).
- :bug: Fix cancer report bugs ([pr48](https://github.com/umccr/gpgr/pull/48)).
- :star: Added CLI support for umccrise cancer report.
  - Modularised CLI subparsers into separate files
- :wrench: Refactored cancer report to use `{sigrap}` pkg for signature tools.
- :wrench: factorise `ClusterId` in LINX VisSvData table ([pr46](https://github.com/umccr/gpgr/pull/46)).
- :wrench: subset multiqc columns ([commit-b7331db](https://github.com/umccr/gpgr/commit/b7331db86151bcbd82268dad264acf823df7471d)).

# gpgr 1.2.0 (2022-02-07)

- :star: Support for DVC ([pr43](https://github.com/umccr/gpgr/pull/43), [issue41](https://github.com/umccr/gpgr/issues/41)).
- :wrench: Fixes for LINX tables ([pr40](https://github.com/umccr/gpgr/pull/40)).
- :wrench: Bug fix for LINX CLI ([8e89ac6](https://github.com/umccr/gpgr/commit/8e89ac67ba45d64e814772b1c12d6fc3b8e7a45d)).

# gpgr 1.1.0 (2022-01-12)

- :star: Add LINX RMarkdown report.
- :star: Added script for tidying umccrised MultiQC JSON output ([pr38](https://github.com/umccr/gpgr/pull/38)).
- :star: Add CLI for LINX reporter ([pr39](https://github.com/umccr/gpgr/pull/39)).
- :wrench: Add LINX description tibble and refactor LINX code ([pr36](https://github.com/umccr/gpgr/pull/36)).

# gpgr 1.0.0 (2021-12-22)

- **Major change**: Moving signature tool wrappers (CHORD, HRDetect, and
  MutationalPatterns) to sigrap (https://github.com/umccr/sigrap)
  ([pr35](https://github.com/umccr/gpgr/pull/35)).
- Support for LINX tables ([pr34](https://github.com/umccr/gpgr/pull/34)).

# gpgr 0.2.0 (2021-11-23)

- new release after a year in the freezer
- :star: add CLI for CHORD, HRDetect and MutationalPatterns
  ([pr26](https://github.com/umccr/gpgr/pull/26),
  [pr28](https://github.com/umccr/gpgr/pull/28),
  [issue17](https://github.com/umccr/gpgr/issues/17))
- :star: add 6 new SBS signatures from the March 2021 COSMIC release (v3.2)
  ([issue22](https://github.com/umccr/gpgr/issues/22),
  [pr23](https://github.com/umccr/gpgr/pull/23))
- :star: update to use the base pipe operator (`|>`) introduced in R 4.1
  ([pr24](https://github.com/umccr/gpgr/pull/24))
- :star: add handy `get_genome_obj` function that returns a BSgenome object, provided
  a 'hgXX' string (useful for MutationalPatterns & CHORD)
  ([pr20](https://github.com/umccr/gpgr/pull/20))
- :star: add handy `pkg_exists` function that checks if the specified R package
  exists on the local system ([pr20](https://github.com/umccr/gpgr/pull/20))
- :question: add `readr::local_edition(1)` before `readr::read_tsv` due to
  `readr` using `vroom` under the hood (I'm not too comfortable with this,
  but anyway). See <https://www.tidyverse.org/blog/2021/07/readr-2-0-0/#readr-2nd-edition>.
  I think it happens because I'm using the `comment` parameter to read in the
  gzipped VCF. Oh well ([pr21](https://github.com/umccr/gpgr/pull/21))
- :wrench: use [bumpversion](https://github.com/c4urself/bump2version) for
  updating pkg versions
  ([issue18](https://github.com/umccr/gpgr/issues/18),
  [pr19](https://github.com/umccr/gpgr/pull/19))
- :wrench: add pre-commit hooks and styler
  ([issue25](https://github.com/umccr/gpgr/issues/25),
  [pr27](https://github.com/umccr/gpgr/pull/27))

# gpgr 0.1.0 (2020-11-19)

- Handle umccrise cancer report functionality:
  - PURPLE CNVs
  - SVs (tables + BND plots)
  - AF (summary + plot)
  - HRD
  - write final tables to tsv + json
  - Support for MutationalPatterns v3.0 (https://github.com/UMCUGenetics/MutationalPatterns)
  - Support for PURPLE v2.51+

# gpgr 0.0.6 (2020-10-09)

- Initial release of gpgr. Features include:
  - Support for CHORD (https://github.com/UMCUGenetics/CHORD)
  - Support for HRDetect (https://github.com/Nik-Zainal-Group/signature.tools.lib)
  - Support for PURPLE (https://github.com/hartwigmedical/hmftools/tree/master/purple)
  - R package webpage deployed via GitHub Actions at https://umccr.github.io/gpgr/
