#!/usr/bin/env bash

set -euo pipefail

gpgr_dir="/Users/pdiakumis/projects/gpgr"
sbj="SBJ02299"
sample="SBJ02299_PRJ221222_L2200695"
dd="${gpgr_dir}/nogit/linx/${sbj}"
pd="${dd}/plots"
td="${dd}/annotations"
od="${gpgr_dir}/nogit/linx/html_reports"
of="${sample}_linx.html"

# link to cli
gpgr_cli=$(Rscript -e 'x = system.file("cli/gpgr.R", package = "gpgr"); cat(x, "\n")')

$gpgr_cli linx \
    --sample ${sample} \
    --plot ${pd} \
    --table ${td} \
    --out ${od}/${of}
