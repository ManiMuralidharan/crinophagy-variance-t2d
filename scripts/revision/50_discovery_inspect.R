# =====================================================================
# 50_discovery_inspect.R
# Prints the structure of the discovery files. Computes nothing.
# RUN:  source("50_discovery_inspect.R")
# =====================================================================

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))

show_tab <- function(path, n = 3) {
  if (!file.exists(path)) { cat("\n== MISSING:", path, "\n"); return(invisible()) }
  d <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  cat("\n==", path, "  (", nrow(d), "rows x", ncol(d), "cols )\n")
  cat("columns:\n"); print(names(d))
  cat("first", n, "rows:\n"); print(head(d, n), row.names = FALSE)
  # flag the columns that matter
  hits <- grep("group|diagnos|glycem|status|t2d|batch|age|sex|bmi|imm|beta|donor|id",
               names(d), ignore.case = TRUE, value = TRUE)
  if (length(hits)) { cat("candidate key columns:\n"); print(hits) }
  invisible(d)
}

show_tab("modelling_dataset_py.csv")
show_tab("donor.csv")
show_tab("ssgsea_scores_with_metadata.csv")

for (f in c("proteomics/proc_prot.rds", "proteomics/unproc_prot.rds")) {
  cat("\n==", f, "\n")
  if (!file.exists(f)) { cat("MISSING\n"); next }
  x <- readRDS(f)
  cat("class:", paste(class(x), collapse = "/"), "\n")
  if (is.matrix(x) || is.data.frame(x)) {
    cat("dim:", paste(dim(x), collapse = " x "), "\n")
    cat("rownames[1:5]:", paste(head(rownames(x), 5), collapse = ", "), "\n")
    cat("colnames[1:5]:", paste(head(colnames(x), 5), collapse = ", "), "\n")
    cat("any NA:", anyNA(x), "  value range:",
        paste(round(range(x, na.rm = TRUE), 2), collapse = " to "), "\n")
  } else {
    str(x, max.level = 2)
  }
  rm(x); gc(verbose = FALSE)
}

# Is GABARAP present in the discovery proteomics, under either key?
cat("\n== GABARAP lookup ==\n")
for (f in c("proteomics/proc_prot.rds", "proteomics/unproc_prot.rds")) {
  if (!file.exists(f)) next
  x <- readRDS(f)
  rn <- rownames(x)
  cat(f, ": rowname style =", paste(head(rn, 2), collapse = ", "), "\n")
  # Entrez for GABARAP is 11337; UniProt is O95166
  for (k in c("11337", "O95166", "GABARAP"))
    cat("   ", k, ":", k %in% rn, "\n")
  rm(x); gc(verbose = FALSE)
}

cat("\nAlso paste: read.csv('C1_C2_between_set.csv')\n")
