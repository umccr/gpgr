## code to prepare `EFFECT_ABBREVIATIONS` dataset goes here

EFFECT_ABBREVIATIONS <- c(
  "3_prime_UTR_truncation" = "3UTRtrunc", "3_prime_UTR_variant" = "3UTRvar",
  "5_prime_UTR_truncation" = "5UTRtrunc",  "5_prime_UTR_variant" = "5UTRvar",
  "feature_fusion" = "Fus", "bidirectional_gene_fusion" = "BidFusG", "gene_fusion" = "FusG",
  "chromosome_number_variation" = "ChromNumV", "conservative_inframe_deletion" = "ConsInframeDel",
  "downstream_gene_variant" = "DnstreamGV", "upstream_gene_variant" = "UpstreamGV",
  "duplication" = "Dup", "exon_loss_variant" = "ExonLossV",
  "feature_ablation" = "DelG", "transcript_ablation" = "DelTx",
  "frameshift_variant" = "FrameshiftV", "intergenic_region" = "IntergenReg", "intragenic_variant" = "IntragenV",
  "intron_variant" = "IntronV",  "no_func_effect" = "NoFuncEff", "no_prio_effect" = "NoPrioEff",
  "non_coding_transcript_variant" = "NoncodTxV",
  "splice_acceptor_variant" = "SpliceAccV",
  "splice_donor_variant" = "SpliceDonV", "splice_region_variant" = "SpliceRegV",
  "start_lost" = "StartLoss", "stop_gained" = "StopGain", "stop_lost" = "StopLoss",
  "TF_binding_site_variant" = "TFBSVar",  "TFBS_ablation" = "TFBSDel")

usethis::use_data(EFFECT_ABBREVIATIONS, overwrite = TRUE)
