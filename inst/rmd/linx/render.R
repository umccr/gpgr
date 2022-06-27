#!/usr/bin/env Rscript

require(here)
require(gpgr)
require(glue)

sbj <- "SBJ02299"
sample <- "SBJ02299_PRJ221222_L2200695"
dd <- glue(here("nogit/linx/{sbj}"))
pd <- glue("{dd}/plots")
td <- glue("{dd}/annotations")
od <- here("nogit/linx/html_reports")
of <- glue("{sample}_linx.html")

gpgr::linx_rmd(
  sample = sample,
  table_dir = td,
  plot_dir = pd,
  out_file = file.path(od, of)
)
