library(dplyr, include.only = c("bind_rows", "left_join"))
library(purrr, include.only = c("map", "map_chr", "reduce", "list_merge", "imap", "set_names"))
library(RJSONIO, include.only = "fromJSON")

mj2df <- function(json) {
  stopifnot(file.exists(json))
  p <- fromJSON(json)
  nm <- "umccr_subj_id"
  gen <- parse_gen(p) |>
    remove_control_samples() |>
    bind_rows(.id = nm)
  raw <- parse_raw(p) |>
    remove_control_samples() |>
    bind_rows(.id = nm)

  # data is in tidy format:
  # - each variable has its own column
  # - each observation (sample) has its own row
  # - each value has its own cell
  left_join(gen, raw, by = nm, suffix = c(".gen", ".raw"))
}

remove_control_samples <- function(l) {
  controls <- c("Alice", "Bob", "Chen", "Elon", "Dakota")
  controls2 <- paste0(controls, rep(c("_T", "_B"), each = length(controls)))
  # also remove 'idxstats' samples
  controls_regex <- paste(c(controls2, "idxstats"), collapse = "|")
  l <- l[!grepl(controls_regex, names(l))]
  l
}

mkdir <- function(d) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
  }
  d
}

# Following are modified/simplified from TidyMultiqc
parse_gen <- function(p) {
  el <- "report_general_stats_data"
  stopifnot(inherits(p, "list"), el %in% names(p))
  p[[el]] |> reduce(~ list_merge(.x, !!!.y))
}

parse_raw <- function(p) {
  el <- "report_saved_raw_data"
  stopifnot(inherits(p, "list"), el %in% names(p))
  # For each tool
  p[[el]][["dragen_frag_len"]] <- NULL
  p[[el]] |>
    imap(function(samples, tool) {
      # Remove the superflous "multiqc_" from the start of the tool name
      tool <- sub("multiqc_", "", tool)

      # For each sample
      samples |> kv_map(function(metrics, sample) {
        # For each metric in the above tool
        list(
          key = sample,
          value = metrics |> kv_map(function(mvalue, mname) {
            # Sanitise metric names
            combined_metric <- list(
              key = paste0(tool, ".", mname),
              value = mvalue
            )
          })
        )
      })
    }) |>
    reduce(utils::modifyList)
}

kv_map <- function(l, func) {
  mapped <- imap(l, func) |>
    set_names(nm = NULL)
  keys <- mapped |> map_chr("key")
  vals <- mapped |> map("value")
  vals |> set_names(keys)
}

select_column_subset_alignmentqc <- function(d) {
  cols_to_keep <- c(
    "umccr_subj_id"           = "umccr_subj_id",
    "tot_input_reads"         = "dragen_map_metrics.Total input reads",
    "tot_input_reads_pct"     = "dragen_map_metrics.Total input reads pct",
    "num_dup_reads"           = "dragen_map_metrics.Number of duplicate marked reads",
    "num_dup_reads_pct"       = "dragen_map_metrics.Number of duplicate marked reads pct",
    "num_uniq_reads_pct"      = "dragen_map_metrics.Number of unique reads (excl. duplicate marked reads) pct",
    "mapped_reads"            = "dragen_map_metrics.Mapped reads",
    "mapped_reads_pct"        = "dragen_map_metrics.Mapped reads pct",
    "unmapped_reads"          = "dragen_map_metrics.Unmapped reads",
    "unmapped_reads_pct"      = "dragen_map_metrics.Unmapped reads pct",
    "singleton_reads"         = "dragen_map_metrics.Singleton reads (itself mapped; mate unmapped)",
    "singleton_reads_pct"     = "dragen_map_metrics.Singleton reads (itself mapped; mate unmapped) pct",
    "paired_reads"            = "dragen_map_metrics.Paired reads (itself & mate mapped)",
    "paired_reads_pct"        = "dragen_map_metrics.Paired reads (itself & mate mapped) pct",
    "paired_reads_proper"     = "dragen_map_metrics.Properly paired reads",
    "paired_reads_proper_pct" = "dragen_map_metrics.Properly paired reads pct",
    "read_length"             = "dragen_map_metrics.Estimated read length",
    "insert_length_mean"      = "dragen_map_metrics.Insert length: mean",
    "insert_length_median"    = "dragen_map_metrics.Insert length: median",
    "contamination"           = "dragen_map_metrics.Estimated sample contamination",
    "cov_seqed_avg_genome"    = "dragen_map_metrics.Average sequenced coverage over genome",
    "cov_aligned_avg_genome"  = "dragen_cov_metrics.Average alignment coverage over genome",
    "cov_autosomal_median"    = "dragen_ploidy.Autosomal median coverage",
    "cov_chrx_median"         = "dragen_ploidy.X median coverage",
    "cov_chry_median"         = "dragen_ploidy.Y median coverage",
    "cov50_genome_pct"        = "dragen_cov_metrics.PCT of genome with coverage [ 50x: inf)",
    "ploidy_est"              = "dragen_ploidy.Ploidy estimation"
  )

  assertthat::assert_that(base::all(cols_to_keep %in% base::names(d)))

  d |>
    dplyr::select(dplyr::all_of(cols_to_keep)) |>
    purrr::set_names(base::names(cols_to_keep))
}
