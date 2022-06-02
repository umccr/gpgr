linx_add_args <- function(subp) {
  linx <- subp$add_parser("linx", help = "UMCCR LINX Report.")
  linx$add_argument("--sample", help = "Sample name.", required = TRUE)
  linx$add_argument("--plot", help = "Path to LINX plot directory.", required = TRUE)
  linx$add_argument("--table", help = "Path to LINX table directory.", required = TRUE)
  linx$add_argument("--out", help = "HTML output file name [def: linx_{sample}.html].")
  linx$add_argument("--quiet", help = "Suppress log printing during rendering.", action = "store_true")
}

linx_parse_args <- function(args) {
  # print(c("You've called linx. Here are the arguments: ", args))
  cli::cli_h1("Start rendering LINX R Markdown report!")
  tab <- normalizePath(args$table)
  pl <- normalizePath(args$plot)
  res <- gpgr::linx_rmd(
    sample = args$sample,
    table_dir = tab,
    plot_dir = pl,
    out_file = args$out,
    quiet = args$quiet
  )
  cli::cli_h1("Finished rendering LINX R Markdown report!")
  cli::cli_alert_info("Path to HTML output:\n{res}")
}
