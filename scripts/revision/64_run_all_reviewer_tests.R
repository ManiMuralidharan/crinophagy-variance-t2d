# Run all four IJMS reviewer validation scripts from RStudio.
# Place this file alongside scripts 60-63 in the Human Islets folder.
# Each script is isolated in its own environment; errors do not prevent
# the remaining scripts from running. Do not treat a partial run as complete.

ROOT <- "."
if (!dir.exists(ROOT)) stop("Project folder not found: ", ROOT)
setwd(ROOT)
OUT <- file.path(ROOT, "reviewer_validation_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

Sys.setenv(REVIEWER_NBOOT = "1200", REVIEWER_NSIM = "350",
           REVIEWER_COVER_BOOT = "300")

jobs <- c("60_reviewer_replication_validation.R",
          "62_discovery_coverage_audit.R",
          "63_bootstrap_coverage_sensitivity.R",
          "61_paired_omics_check.R")

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("Missing 'readxl'. Run install.packages('readxl') first.")
}
if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
    !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  warning("Bioconductor mapping packages not installed. Discovery mapping in script 62 may be incomplete. Install with BiocManager::install(c('AnnotationDbi','org.Hs.eg.db')).")
}

status <- data.frame(script = jobs, status = NA_character_, detail = NA_character_,
                     stringsAsFactors = FALSE)

for (i in seq_along(jobs)) {
  f <- jobs[i]
  cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat("RUNNING ", f, "\n", sep = "")
  cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")

  if (!file.exists(f)) {
    status$status[i] <- "MISSING SCRIPT"
    status$detail[i] <- paste("Copy", f, "into", ROOT)
    next
  }
  if (identical(f, "61_paired_omics_check.R") &&
      !file.exists(file.path(ROOT, "paired_omics_long.csv"))) {
    # Execute once to create the expected input template, then flag as pending.
    tryCatch(source(f, local = new.env(parent = globalenv())),
             error = function(e) cat("Paired-data template created (if possible): ",
                                    conditionMessage(e), "\n", sep = ""))
    status$status[i] <- "INPUT NEEDED"
    status$detail[i] <- "Prepare paired_omics_long.csv from the input template, then rerun script 61"
    next
  }

  outcome <- tryCatch({
    source(f, local = new.env(parent = globalenv()))
    list(ok = TRUE, message = "Completed; inspect the output CSVs and console warnings")
  }, error = function(e) {
    list(ok = FALSE, message = conditionMessage(e))
  })
  status$status[i] <- if (outcome$ok) "RAN" else "FAILED"
  status$detail[i] <- outcome$message
  if (outcome$ok && identical(f, "62_discovery_coverage_audit.R")) {
    raw_file <- file.path(OUT, "10_unprocessed_coverage_status.txt")
    if (file.exists(raw_file)) {
      raw_state <- readLines(raw_file, n=1, warn=FALSE)
      if (length(raw_state) && !identical(raw_state, "unprocessed_matrix_parsed")) {
        status$status[i] <- "PARTIAL - RAW DATA"
        status$detail[i] <- paste("Processed coverage audited; raw detection unresolved:", raw_state)
      }
    }
  }
  cat("STATUS: ", status$status[i], " | ", status$detail[i], "\n", sep = "")
}

file <- file.path(OUT, "00_run_status.csv")
write.csv(status, file, row.names = FALSE)
cat("\n================== FINAL STATUS ==================\n")
print(status, row.names = FALSE, right = FALSE)
cat("\nSaved status: ", file, "\n", sep = "")
cat("Results folder: ", OUT, "\n", sep = "")
cat("'RAN' means only that the script finished, not that the results have been scientifically validated.\n")

# The 60/62 scripts audit processed data even if raw coverage cannot be parsed.
# Check 10_unprocessed_coverage_status.txt before claiming raw detection rates.
