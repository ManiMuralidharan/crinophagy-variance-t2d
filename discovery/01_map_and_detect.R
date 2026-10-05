# =====================================================================
# Step 1 (corrected): Entrez mapping + detection frequency
#
# Findings from step 0 that this script acts on:
#   - Feature IDs are NCBI Entrez Gene IDs, not symbols
#   - proc_prot.rds has 0 NAs => imputed or complete-case filtered
#   - unproc_prot.rds has 7.1% NA => this is where missingness lives
#   - 129 proteomics donors, 87 also have bulk RNA
# =====================================================================

# ---- 0. Setup -------------------------------------------------------

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db")
# install.packages(c("dplyr", "ggplot2"))

library(org.Hs.eg.db)
library(AnnotationDbi)
library(dplyr)
library(ggplot2)

setwd("/Users/mmani3/Downloads/Human Islets")

prot_proc <- readRDS("./proteomics/proc_prot.rds")
prot_raw  <- readRDS("./proteomics/unproc_prot.rds")


# ---- 1. Gene sets (frozen before analysis) --------------------------

crinophagy <- c(
  "UNC13D", "STX7", "STX8", "STX17", "VAMP7", "VAMP8", "SNAP23",
  "VPS33A", "VPS16", "VPS11",
  "RAB3A", "RAB3GAP1", "RAB27A", "CADPS", "SYTL4", "PCLO",
  "LAMP1", "LAMP2", "CTSB", "CTSD", "CTSL",
  "ATP6V0D1", "ATP6V1A", "TFEB", "PIKFYVE",
  "INS", "CHGA", "CHGB", "PCSK1", "PCSK2", "CPE"
)

macroautophagy <- c(
  "MAP1LC3B", "SQSTM1", "ATG5", "ATG7", "BECN1", "WIPI2", "ULK1"
)

core_three <- c("UNC13D", "CADPS", "RAB3A")


# ---- 2. Map symbols -> Entrez ---------------------------------------
# Mapped programmatically rather than hardcoded, so the IDs are correct
# by construction and the version is recorded for the methods section.

map_symbols <- function(symbols) {
  m <- AnnotationDbi::select(org.Hs.eg.db,
                             keys = symbols,
                             keytype = "SYMBOL",
                             columns = "ENTREZID")
  m[!is.na(m$ENTREZID), ]
}

map_crino <- map_symbols(crinophagy)
map_macro <- map_symbols(macroautophagy)

cat("===== SYMBOL -> ENTREZ MAPPING =====\n")
print(map_crino)

unmapped <- setdiff(crinophagy, map_crino$SYMBOL)
if (length(unmapped)) {
  cat("\nUnmapped symbols (check for aliases):",
      paste(unmapped, collapse = ", "), "\n")
}
cat("\norg.Hs.eg.db version:",
    as.character(packageVersion("org.Hs.eg.db")), "\n\n")


# ---- 3. Did the targets survive processing? -------------------------
# This is the real feasibility question. 2056 proteins were dropped
# during filtering, most likely the high-dropout ones.

ids_proc <- as.character(prot_proc[[1]])
ids_raw  <- as.character(prot_raw[[1]])

survival <- map_crino %>%
  mutate(
    in_raw  = ENTREZID %in% ids_raw,
    in_proc = ENTREZID %in% ids_proc,
    fate = case_when(
      !in_raw            ~ "never measured",
      in_raw & !in_proc  ~ "measured, filtered out",
      TRUE               ~ "retained"
    )
  )

cat("===== TARGET PROTEIN FATE =====\n")
print(survival[, c("SYMBOL", "ENTREZID", "fate")])

cat("\nSummary:\n")
print(table(survival$fate))
cat("\n")


# ---- 4. Detection frequency in the RAW matrix -----------------------

raw_mat <- as.matrix(prot_raw[, -1, drop = FALSE])
storage.mode(raw_mat) <- "numeric"
rownames(raw_mat) <- ids_raw

n_donors <- ncol(raw_mat)

# Detection function supporting either SYMBOL or ENTREZID matching
detection <- function(symbols_vec) {
  # Match against rownames (works whether rownames are symbols or Entrez IDs)
  present_ids <- intersect(symbols_vec, rownames(raw_mat))
  
  # Fallback: if rownames are Entrez IDs, try matching mapped Entrez IDs
  if (length(present_ids) == 0) {
    mapped_entrez <- map_crino$ENTREZID[map_crino$SYMBOL %in% symbols_vec]
    present_ids <- intersect(mapped_entrez, rownames(raw_mat))
  }
  
  if (length(present_ids) == 0) return(NULL)
  
  data.frame(
    id = present_ids,
    n_detected = sapply(present_ids, function(id) {
      v <- raw_mat[id, ]
      sum(!is.na(v) & v != 0)
    }),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      symbol = ifelse(id %in% map_crino$SYMBOL, id, 
                      map_crino$SYMBOL[match(id, map_crino$ENTREZID)]),
      pct_detected = round(100 * n_detected / n_donors, 1)
    ) %>%
    arrange(desc(pct_detected))
}

det_crino <- detection(crinophagy)
det_macro <- detection(macroautophagy)

cat("===== DETECTION FREQUENCY (raw matrix, n =", n_donors, ") =====\n")
print(det_crino)
cat("\n===== COMPARATOR SET =====\n")
print(det_macro)


# ---- 5. The decision ------------------------------------------------

if (is.null(det_crino) || nrow(det_crino) == 0) {
  cat("\nNone of the target proteins were found/detected in the raw matrix.\n")
  cat("=> SWITCH to transcript-primary design.\n")
} else {
  core <- det_crino[det_crino$symbol %in% core_three, ]
  
  if (is.null(core) || nrow(core) == 0) {
    cat("\nNone of UNC13D / CADPS / RAB3A measured in core machinery.\n")
    cat("=> SWITCH to transcript-primary design.\n")
  } else {
    min_core <- min(core$pct_detected)
    cat("\nLowest core detection:", min_core, "%\n")
    if (min_core >= 70) {
      cat("=> PROCEED with protein-primary design.\n")
    } else if (min_core >= 50) {
      cat("=> BORDERLINE. Inspect effect sizes before committing.\n")
    } else {
      cat("=> SWITCH to transcript-primary, protein as confirmation.\n")
    }
  }
}


# ---- 6. Save the record ---------------------------------------------

write.csv(survival,  "target_protein_fate.csv",     row.names = FALSE)
write.csv(det_crino, "detection_crinophagy.csv",    row.names = FALSE)
write.csv(det_macro, "detection_macroautophagy.csv", row.names = FALSE)

if (!is.null(det_crino)) {
  det_crino %>%
    mutate(symbol = factor(symbol, levels = rev(symbol))) %>%
    ggplot(aes(x = pct_detected, y = symbol)) +
    geom_col(fill = "#2e7d32") +
    geom_vline(xintercept = c(50, 70), linetype = "dashed") +
    labs(x = "Donors with quantified value (%)", y = NULL,
         title = "Crinophagy protein detection",
         subtitle = paste("HumanIslets raw proteomics, n =", n_donors)) +
    theme_minimal(base_size = 11)
  ggsave("detection_crinophagy.pdf", plot = p, width = 6, height = 8)
}
