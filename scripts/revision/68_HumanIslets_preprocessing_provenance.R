# 68_HumanIslets_preprocessing_provenance.R
# IJMS Referee 2: provenance and export-matrix reconciliation audit.
# This does NOT claim to reproduce the provider's original missForest fit.
# Requires only base R and 15_raw_uniprot_entrez_mapping.csv from script 65.
# Run: source("68_HumanIslets_preprocessing_provenance.R")
# Files are read only; outputs go to reviewer_validation_output.

ROOT <- "."
OUT <- file.path(ROOT, "reviewer_validation_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

paths <- c(
  raw = file.path(ROOT, "proteomics", "unproc_prot.rds"),
  processed = file.path(ROOT, "proteomics", "proc_prot.rds"),
  donor = file.path(ROOT, "donor.csv"),
  mapping = file.path(OUT, "15_raw_uniprot_entrez_mapping.csv")
)
missing_files <- paths[!file.exists(paths)]
if (length(missing_files)) stop("Missing inputs: ", paste(missing_files, collapse = " | "))

read_feature_matrix <- function(path, id_col) {
  d <- readRDS(path)
  if (!is.data.frame(d) || !id_col %in% names(d))
    stop("Unexpected data format in ", path, "; expected id column ", id_col)
  ids <- as.character(d[[id_col]])
  if (anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids)))
    stop("Missing or duplicate feature IDs in ", path)
  m <- as.matrix(d[, setdiff(names(d), id_col), drop = FALSE])
  storage.mode(m) <- "numeric"
  rownames(m) <- ids
  if (anyDuplicated(colnames(m))) stop("Duplicate donor columns in ", path)
  m
}

RAW <- read_feature_matrix(paths[["raw"]], "feature_id")
PRO <- read_feature_matrix(paths[["processed"]], "gene_id")
D <- read.csv(paths[["donor"]], check.names = FALSE, stringsAsFactors = FALSE)
mp <- read.csv(paths[["mapping"]], check.names = FALSE, stringsAsFactors = FALSE,
               colClasses = "character")
if (!all(c("feature_id", "mapped_entrez", "status") %in% names(mp)))
  stop("Mapping file must contain feature_id, mapped_entrez and status")
if (anyDuplicated(mp$feature_id)) stop("Duplicate accessions in mapping file")
if (!setequal(mp$feature_id, rownames(RAW)))
  stop("Mapping-file accessions differ from raw matrix; regenerate script 65 results")
if (!all(c("record_id", "diagnosis_computed") %in% names(D)))
  stop("Donor metadata missing identifier or diagnosis fields")
if (anyDuplicated(D$record_id)) stop("Duplicate donor records")

# Data provenance: cryptographic file identity can be checked by a collaborator;
# filesystem timestamps are NOT a reliable original download date.
manifest <- do.call(rbind, lapply(names(paths), function(k) {
  p <- paths[[k]]
  fi <- file.info(p)
  data.frame(role = k, path = normalizePath(p), bytes = as.numeric(fi$size),
             filesystem_modified = format(fi$mtime, "%Y-%m-%d %H:%M:%S %Z"),
             md5 = unname(tools::md5sum(p)), stringsAsFactors = FALSE)
}))
write.csv(manifest, file.path(OUT, "30_export_file_manifest.csv"), row.names = FALSE)

shared <- intersect(colnames(RAW), colnames(PRO))
if (length(shared) != 129L)
  stop("Expected 129 shared donors; got ", length(shared), ". Stop to verify export version")
if (!all(shared %in% D$record_id)) stop("Some omics donor IDs absent from donor.csv")

D$grp <- unname(c(None = "ND", `Pre.T2D` = "preT2D", Type2 = "T2D")[
  trimws(as.character(D$diagnosis_computed))])
COV <- c("donorage", "donorsex", "bodymassindex",
         "predistributionculturetime", "percentieqrecoverypostculture")
if (!all(COV %in% names(D)))
  stop("Discovery covariates missing: ", paste(setdiff(COV, names(D)), collapse = ", "))
D <- D[D$record_id %in% shared & D$grp %in% c("ND", "preT2D", "T2D") &
       !is.na(D$grp) & complete.cases(D[, COV, drop = FALSE]), , drop = FALSE]
D$donorsex <- factor(D$donorsex)
D$grp <- factor(D$grp, levels = c("ND", "preT2D", "T2D"))
expected <- c(ND = 87L, preT2D = 20L, T2D = 16L)
actual <- table(D$grp)
print(actual)
if (!identical(as.integer(actual), as.integer(expected)))
  stop("Full covariate-complete groups do not reproduce expected 87/20/16")

# 50% feature-missingness filter is assessed on ALL 129 donor columns.
# This is a structural comparison only: exact source stage order and dynamic
# cross-release sample lists are not established by the present exports.
raw_missing_n <- rowSums(!is.finite(RAW[, shared, drop = FALSE]))
raw_missing_prop <- raw_missing_n / length(shared)
raw_status <- mp$status[match(rownames(RAW), mp$feature_id)]
raw_entrez <- mp$mapped_entrez[match(rownames(RAW), mp$feature_id)]
matched_entrez <- !is.na(raw_entrez) & nzchar(raw_entrez) & raw_status == "UNIQUE"
flag_in_processed <- matched_entrez & raw_entrez %in% rownames(PRO)
feature_audit <- data.frame(
  feature_id = rownames(RAW), entrez = raw_entrez,
  annotation_status = raw_status,
  n_donors_with_measurement = length(shared) - raw_missing_n,
  missing_fraction_129 = raw_missing_prop,
  pass_documented_at_most_50pct_missing = raw_missing_prop <= .50,
  maps_uniquely_to_entrez = matched_entrez,
  mapped_gene_retained_in_processed = flag_in_processed,
  stringsAsFactors = FALSE
)
write.csv(feature_audit, file.path(OUT, "31_raw_to_processed_structural_audit.csv"), row.names = FALSE)
cat("\nRaw features:", nrow(RAW), " Processed features:", nrow(PRO), "\n")
cat("Raw features passing <=50% missingness:", sum(raw_missing_prop <= .50), "\n")
cat("Unique Entrez mapped raw features:", sum(matched_entrez), "\n")
cat("Mapped raw entries whose Entrez is in the processed matrix:",
    sum(flag_in_processed), "\n")
cat("Raw features with >50% missing values but gene retained processed:",
    sum(raw_missing_prop > .50 & flag_in_processed), "\n")

# UNC13D must have precisely one unambiguous raw mapping in the present audit.
uni_unc <- mp$feature_id[mp$mapped_entrez == "201294" & mp$status == "UNIQUE"]
uni_unc <- intersect(uni_unc, rownames(RAW))
if (length(uni_unc) != 1L || !"201294" %in% rownames(PRO))
  stop("UNC13D accession/Entrez is ambiguous or not in processed data")
cat("UNC13D raw accession:", uni_unc, "\n")

# Compare several model fits to isolate a key Methods distinction:
# manuscript states 123 donors used for nuisance adjustment and 103 donors
# used for the ND-vs-T2D SD-ratio computation.
# The previously provided script 67 instead refitted on 103 or 74 donors.
# This comparison is descriptive, not a fully calibrated variance test.
all_ids <- colnames(RAW)
raw_intensity <- RAW[uni_unc, all_ids]
raw_intensity[!is.finite(raw_intensity) | raw_intensity <= 0] <- NA_real_
raw_for_norm <- RAW[, all_ids, drop = FALSE]
raw_for_norm[!is.finite(raw_for_norm) | raw_for_norm <= 0] <- NA_real_
raw_log <- log2(raw_for_norm)
sample_log_median <- apply(raw_log, 2, median, na.rm = TRUE)
if (any(!is.finite(sample_log_median))) stop("Invalid donor normalization medians")
# Median-centering plus a common constant is a sensitivity approximation only.
# It has not been verified to exactly match HumanIslets provider processing.
raw_unc_log <- log2(raw_intensity)
raw_unc_centered <- raw_unc_log - sample_log_median + median(sample_log_median)
proc_unc <- PRO["201294", all_ids]

estimate <- function(label, values, fit_selection = c("all123", "two103", "raw74")) {
  fit_selection <- match.arg(fit_selection)
  d <- D
  d$y <- as.numeric(values[as.character(d$record_id)])
  if (fit_selection == "two103") d <- d[d$grp %in% c("ND", "T2D"), , drop = FALSE]
  if (fit_selection == "raw74") {
    d <- d[d$grp %in% c("ND", "T2D") & is.finite(d$y), , drop = FALSE]
  }
  d <- d[is.finite(d$y), , drop = FALSE]
  d$donorsex <- droplevels(factor(d$donorsex))
  X <- model.matrix(~ donorage + donorsex + bodymassindex +
                     predistributionculturetime + percentieqrecoverypostculture, data = d)
  if (qr(X)$rank < ncol(X)) stop("Rank deficient design in ", label)
  fit <- lm(y ~ donorage + donorsex + bodymassindex +
              predistributionculturetime + percentieqrecoverypostculture, data = d)
  rr <- residuals(fit)
  nd <- d$grp == "ND"; td <- d$grp == "T2D"
  if (sum(nd) < 5 || sum(td) < 5) stop("Too few ND/T2D in ", label)
  data.frame(analysis = label, fit_population = fit_selection,
             n_fit = nrow(d), n_ND_ratio = sum(nd), n_T2D_ratio = sum(td),
             sd_ND = sd(rr[nd]), sd_T2D = sd(rr[td]),
             sd_ratio = sd(rr[td]) / sd(rr[nd]),
             stringsAsFactors = FALSE)
}

results <- rbind(
  estimate("processed_full_fit123", proc_unc, "all123"),
  estimate("processed_fit103", proc_unc, "two103"),
  estimate("processed_fit74_raw_observed", replace(proc_unc, !is.finite(raw_unc_log), NA_real_), "raw74"),
  estimate("raw_log2_fit_on_measured_all_groups", raw_unc_log, "all123"),
  estimate("raw_median_centered_fit_on_measured_all_groups", raw_unc_centered, "all123"),
  estimate("raw_log2_fit74", raw_unc_log, "raw74"),
  estimate("raw_median_centered_fit74", raw_unc_centered, "raw74")
)
write.csv(results, file.path(OUT, "32_unc13d_fit_population_reconciliation.csv"), row.names = FALSE)
cat("\nUNC13D different donor fitting populations (descriptive):\n")
print(results, row.names = FALSE, digits = 4)

# How similar are the raw and processed observations? This does NOT
# demonstrate provider pipeline identity or validate imputed estimates.
ids <- as.character(D$record_id)
rawv <- raw_unc_centered[ids]; prov <- proc_unc[ids]
obs <- is.finite(rawv) & is.finite(prov)
qc <- data.frame(record_id = ids, group = as.character(D$grp),
                 raw_observed = is.finite(raw_unc_log[ids]),
                 raw_log2 = unname(raw_unc_log[ids]),
                 raw_median_centered = unname(rawv),
                 processed_value = unname(prov))
write.csv(qc, file.path(OUT, "33_unc13d_raw_processed_donor_comparison.csv"), row.names = FALSE)
obs_r <- if (sum(obs) > 3) cor(rawv[obs], prov[obs]) else NA_real_
cat("UNC13D raw-observed vs processed Pearson r:", round(obs_r, 4),
    "over", sum(obs), "covariate-complete donors including pre-T2D.\n")

# Record visible local documentation without inventing a source-export version.
readme_files <- list.files(file.path(ROOT, "proteomics"),
                           pattern = "readme|manifest|version|processing|metadata",
                           ignore.case = TRUE, full.names = TRUE, recursive = TRUE)
notes <- c(
  "HumanIslets discovery export provenance audit",
  paste("Audit run:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "Published pipeline (HumanIslets paper/technical documentation): >50% missing filter; median normalize; log2; missForest (imp4p); UniProt->Entrez mean aggregation; bottom 20% variance filter.",
  "Provider processing pipeline APPLICABILITY to this exact 129-donor export: NOT VERIFIED.",
  "Original export download timestamp/release/version: NOT VERIFIED by file modification times.",
  "Source technical replicate handling before donor-level export: NOT VERIFIED.",
  "Proteomics MS acquisition batch adjustment: NOT VERIFIED; no batch field was identified in prior donor audit.",
  "Provider's precise median normalization and missForest settings/seed: NOT VERIFIED.",
  "Complete processed values for raw-missing UNC13D donors cannot be labeled proven-imputed without source-lineage confirmation.",
  "Unprocessed export contains donor-level intensity values, not raw mass-spectrometry spectra.",
  "Gene-level differences may arise from filtering, mapping, aggregation, or version changes; numbers alone do not identify a step.",
  "Published discovery UNC13D SD ratio 2.065 was stated as unprocessed per-protein result; do not label processed sensitivity result 2.083 as exact reproduction.",
  "Raw observed 59 ND/15 T2D; processed full 87 ND/16 T2D. Missingness can be informative.",
  "Brown-Forsythe on covariate residuals and percentile bootstrap require calibration; this script computes descriptive SD ratios only.",
  "For exact provenance: obtain timestamped HumanIslets download/analysis bundle; identify matching HumanIsletsR commit and data release; verify source preprocessing and per-protein/technical replicate handling.",
  paste("Local proteomics documentation files:", if (length(readme_files)) paste(readme_files, collapse = " | ") else "NONE")
)
writeLines(notes, file.path(OUT, "34_provenance_interpretation_and_open_questions.txt"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "35_R_sessionInfo.txt"))
cat("\nDONE: 30_manifest, 31_structural_audit, 32_fit_reconciliation,\n",
    "33_donor_comparison, 34_notes, 35_R_sessionInfo in ", OUT, "\n", sep = "")
