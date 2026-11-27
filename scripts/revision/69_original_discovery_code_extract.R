# IJMS historical discovery-analysis code audit (read-only)
# Purpose: extract the lines surrounding early source-code statements that
# actually load protein matrices, fit models, compute SD ratios, and write results.
# Does NOT run the historical analysis and does not modify the source scripts.

ROOT <- "."
OUT  <- file.path(ROOT, "reviewer_validation_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

chosen <- c(
  "00_inspect_download.R",
  "01_map_and_detect.R",
  "02_model_crinophagy.R",
  "03_pipeline.py",
  "04_robustness.py",
  "05_unc13d_check.py",
  "06_rab3a_and_transcript.py",
  "09_replication_gse164416.R",
  "20_layer_concordance.R",
  "21_mean_volcano.R",
  "36_fig5A.R"
)

files <- file.path(ROOT, chosen)
present <- file.exists(files)
if (!any(present)) stop("None of the expected early scripts was found in ", ROOT)

meta <- do.call(rbind, lapply(seq_along(files), function(i) {
  f <- files[i]
  data.frame(
    script = chosen[i], exists = present[i],
    size_bytes = if (present[i]) as.numeric(file.info(f)$size) else NA_real_,
    modified_local = if (present[i]) format(file.info(f)$mtime, "%Y-%m-%d %H:%M:%S %Z") else NA_character_,
    md5 = if (present[i]) unname(tools::md5sum(f)) else NA_character_,
    stringsAsFactors = FALSE
  )
}))
write.csv(meta, file.path(OUT, "37_historical_script_manifest.csv"), row.names = FALSE)

patterns <- list(
  input_and_preprocessing = paste(c(
    "proc_prot", "unproc_prot", "readRDS", "read_rds", "read_r", "pyreadr",
    "feature_id", "gene_id", "Q70J99", "201294", "dropna", "isna",
    "missing", "imput", "log2", "median", "normali", "aggregation"
  ), collapse = "|"),
  model_and_inference = paste(c(
    "sd_ratio", "sd\\(", "std\\(", "levene", "brown", "resid",
    "covar", "lm\\(", "ols\\(", "variance", "fdr", "p.adjust",
    "multipletests", "statsmodels", "benjamini", "nboot", "bootstrap",
    "permut", "2\\.065", "1\\.725"
  ), collapse = "|")
)

report <- file.path(OUT, "37_historical_discovery_code_context.txt")
con <- file(report, open = "wt")
w <- function(...) cat(..., "\n", file = con, sep = "")
w("HISTORICAL SOURCE CODE CONTEXT — READ ONLY")
w("Source folder: ", ROOT)
w("Do not treat file modification dates or comments as evidence of prospective hypothesis registration.")
w("These snippets show code references, not proof that a particular file was run for publication.")

for (f in files[present]) {
  lines <- tryCatch(readLines(f, warn = FALSE), error = function(e) character())
  w("\n", strrep("=", 85))
  w("SCRIPT: ", basename(f), "  (", length(lines), " lines)")
  for (topic in names(patterns)) {
    match_ix <- grep(patterns[[topic]], lines, ignore.case = TRUE, perl = TRUE)
    w("\nTOPIC: ", topic, "  [", length(match_ix), " matching lines]")
    if (!length(match_ix)) next
    if (length(match_ix) > 45L) {
      w("Showing first 45 matches. All line numbers: ", paste(match_ix, collapse = ", "))
      match_ix <- head(match_ix, 45L)
    }
    idx <- sort(unique(unlist(lapply(match_ix, function(k) seq.int(max(1L, k - 3L), min(length(lines), k + 4L))))))
    old <- 0L
    for (j in idx) {
      if (j > old + 1L) w("   ...")
      w(sprintf("%4d  %s", j, lines[j]))
      old <- j
    }
  }
}
close(con)
cat("Wrote:\n", report, "\n",
    file.path(OUT, "37_historical_script_manifest.csv"), "\n")
cat("Open the text report and inspect 02_model_crinophagy.R, 03_pipeline.py, 05_unc13d_check.py, and 06_rab3a_and_transcript.py first.\n")
cat("CAUTION: match context does not establish which script actually generated the published figure or table.\n")
