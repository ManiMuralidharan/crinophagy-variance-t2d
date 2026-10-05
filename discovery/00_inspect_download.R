# =====================================================================
# HumanIslets download: structural inspection
# Run this from the folder that contains the four subfolders.
# Answers: donor n, ID type, overlap, missingness, metadata presence.
# No statistics here. Nothing is tested. This is reconnaissance.
# =====================================================================

setwd("/Users/mmani3/Downloads/Human Islets")   # <- edit this
cat("Working directory:", getwd(), "\n\n")

cat("===== FILES PRESENT =====\n")
print(list.files(recursive = TRUE))
cat("\n")


# ---- Locate the files regardless of exact folder naming -------------

find_file <- function(pattern) {
  hits <- list.files(pattern = pattern, recursive = TRUE, full.names = TRUE)
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

f_prot <- find_file("proc_prot\\.rds$")
f_prot_raw <- find_file("unproc_prot\\.rds$")
f_rna  <- find_file("proc_rnaseq\\.rds$")

cat("proteomics (proc):  ", f_prot, "\n")
cat("proteomics (unproc):", f_prot_raw, "\n")
cat("bulk RNA (proc):    ", f_rna, "\n\n")


# ---- CHECK 1: dimensions and donor count ----------------------------

prot <- readRDS(f_prot)

cat("===== CHECK 1: PROTEOMICS DIMENSIONS =====\n")
cat("Class:", class(prot), "\n")
cat("Dimensions:", dim(prot), "\n\n")

cat("All column names:\n")
print(colnames(prot))
cat("\n")

# Donor IDs are typically numeric-looking codes. Metadata columns are
# words. Eyeball the list above and confirm the split below is right.
meta_cols <- colnames(prot)[1:5]
cat("Assuming first columns are metadata:", meta_cols, "\n")
cat("=> implied donor count:", ncol(prot) - 1, "(if only col 1 is ID)\n\n")


# ---- CHECK 2: what are the feature identifiers? ---------------------

cat("===== CHECK 2: FEATURE ID TYPE =====\n")
cat("First column name:", colnames(prot)[1], "\n")
cat("First 10 values:\n")
print(head(prot[[1]], 10))

# Interpretation:
#   P12345 / Q9NRA8 style  -> UniProt accession, needs mapping
#   INS / LAMP1 / CTSB     -> gene symbol, my gene sets work directly
#   ENSG00000...           -> Ensembl, needs mapping
cat("\n")


# ---- CHECK 3: donor overlap between proteomics and RNA --------------

cat("===== CHECK 3: PROTEIN / RNA DONOR OVERLAP =====\n")

if (!is.na(f_rna)) {
  rna <- readRDS(f_rna)
  cat("RNA dimensions:", dim(rna), "\n")

  # proc_rnaseq has 5 metadata columns per the documentation
  rna_donors  <- colnames(rna)[-(1:5)]
  prot_donors <- colnames(prot)[-1]

  ov <- intersect(prot_donors, rna_donors)

  cat("Donors with proteomics:", length(prot_donors), "\n")
  cat("Donors with bulk RNA:  ", length(rna_donors), "\n")
  cat("OVERLAP (both layers): ", length(ov), "\n")
  cat("This overlap is your real sample size for matched analysis.\n\n")
} else {
  cat("RNA file not found, skipping.\n\n")
}


# ---- CHECK 4: missingness (the imputation question) -----------------

cat("===== CHECK 4: MISSING VALUES =====\n")

num_prot <- prot[, -1, drop = FALSE]
num_prot <- as.matrix(num_prot)
storage.mode(num_prot) <- "numeric"

n_na   <- sum(is.na(num_prot))
n_zero <- sum(num_prot == 0, na.rm = TRUE)
n_tot  <- length(num_prot)

cat("Total cells:  ", n_tot, "\n")
cat("NA values:    ", n_na,   sprintf("(%.1f%%)", 100 * n_na / n_tot), "\n")
cat("Zero values:  ", n_zero, sprintf("(%.1f%%)", 100 * n_zero / n_tot), "\n\n")

if (n_na == 0 && n_zero == 0) {
  cat("WARNING: no NAs and no zeros. This matrix has been imputed.\n")
  cat("Detection frequency is meaningless here. Use unproc_prot.rds\n")
  cat("to assess true missingness, and ask the platform how the\n")
  cat("imputation was done before interpreting any result.\n\n")
} else {
  cat("Missing values present. Detection frequency is interpretable.\n\n")
}

# Compare against the raw matrix
if (!is.na(f_prot_raw)) {
  praw <- readRDS(f_prot_raw)
  nraw <- as.matrix(praw[, -1, drop = FALSE])
  storage.mode(nraw) <- "numeric"
  cat("Raw matrix dims:", dim(praw), "\n")
  cat("Raw NA fraction:",
      sprintf("%.1f%%", 100 * sum(is.na(nraw)) / length(nraw)), "\n")
  cat("Proteins lost in filtering:", nrow(praw) - nrow(prot), "\n\n")
}


# ---- CHECK 5: are the target proteins present at all? ---------------

cat("===== CHECK 5: TARGET PROTEIN PRESENCE =====\n")

core <- c("UNC13D", "CADPS", "RAB3A", "RAB3GAP1", "STX17", "VAMP8",
          "VPS33A", "LAMP1", "LAMP2", "CTSB", "CTSD", "TFEB")

ids <- as.character(prot[[1]])
found <- core[core %in% ids]
missing <- setdiff(core, ids)

cat("Found:  ", if (length(found))   paste(found, collapse = ", ")   else "none", "\n")
cat("Missing:", if (length(missing)) paste(missing, collapse = ", ") else "none", "\n")

if (length(found) == 0) {
  cat("\nNone found by gene symbol. The IDs are probably UniProt or\n")
  cat("Ensembl (see CHECK 2). Mapping needed before this test means\n")
  cat("anything. Do not conclude the proteins are absent.\n")
}
cat("\n")


# ---- CHECK 6: is donor metadata anywhere? ---------------------------

cat("===== CHECK 6: DONOR METADATA =====\n")
csvs <- list.files(pattern = "\\.csv$", recursive = TRUE)
cat("CSV files found:\n"); print(csvs)
cat("\nYou need a table with diabetes status, HbA1c, BMI, age, sex,\n")
cat("culture time and exocrine fraction, keyed by record_id.\n")
cat("If it is not here, download it via 'Download Donor List' on the\n")
cat("Select Donors panel. Without it nothing downstream is analysable.\n")
