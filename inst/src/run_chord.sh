#!/usr/bin/env bash

timestamp=$(date +%s)

./gpgr.R chord \
    --sample "sampleA" \
    --snv "../extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz" \
    --sv "../extdata/umccrise/sv/manta.vcf.gz" \
    --out "../../nogit/chord_results_${timestamp}.json.gz"
