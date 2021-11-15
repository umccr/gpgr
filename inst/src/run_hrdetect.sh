#!/usr/bin/env bash

timestamp=$(date +%s)

./gpgr.R hrdetect \
    --sample "sampleA" \
    --snv "../extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz" \
    --sv "../extdata/umccrise/sv/manta.vcf.gz" \
    --cnv "../extdata/purple/purple.cnv.somatic.tsv" \
    --out "../../nogit/hrdetect_results_${timestamp}.json.gz"
