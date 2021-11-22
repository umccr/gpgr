#!/usr/bin/env bash

timestamp=$(date +%s)

./gpgr.R mutpat \
    --sample "sampleA" \
    --snv "../extdata/umccrise/snv/somatic-ensemble-PASS.vcf.gz" \
    --outdir "../../nogit/mutpat_results"
