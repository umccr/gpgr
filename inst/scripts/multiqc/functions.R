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
  left_join(gen, raw, by = nm)
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
