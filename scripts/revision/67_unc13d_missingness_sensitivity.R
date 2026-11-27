# 67_unc13d_missingness_sensitivity.R
# Referee 2: diagnose UNC13D raw measurement availability and test whether
# discovery SD ratios change when donors lacking raw UNC13D are excluded.
# This is a SENSITIVITY ANALYSIS, not a reconstruction of the original
# proteomics normalization or a validated covariate-adjusted variance test.
# Run AFTER 65_raw_proteomics_uniprot_audit_fixed.R has generated file 15.

ROOT <- "."
OUT  <- file.path(ROOT, "reviewer_validation_output")
BATCH_COL <- ""  # Optional: enter a VERIFIED discovery proteomics batch column
NBOOT <- 1000L
set.seed(20261113)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

read_feature_matrix <- function(path, id_col) {
  if (!file.exists(path)) stop("Missing: ", path)
  x <- readRDS(path)
  if (!is.data.frame(x) || !id_col %in% names(x))
    stop("Expected data.frame with column '", id_col, "' in ", path)
  ids <- as.character(x[[id_col]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids))
    stop("Missing or duplicated feature IDs in ", path)
  x <- as.matrix(x[, setdiff(names(x), id_col), drop = FALSE])
  storage.mode(x) <- "numeric"
  rownames(x) <- ids
  if (anyDuplicated(colnames(x))) stop("Duplicate donor IDs in ", path)
  x
}

RAW <- read_feature_matrix(file.path(ROOT, "proteomics/unproc_prot.rds"), "feature_id")
PRO <- read_feature_matrix(file.path(ROOT, "proteomics/proc_prot.rds"), "gene_id")
map_path <- file.path(OUT, "15_raw_uniprot_entrez_mapping.csv")
if (!file.exists(map_path))
  stop("Missing ", map_path, ". First run 65_raw_proteomics_uniprot_audit_fixed.R")
mp <- read.csv(map_path, stringsAsFactors = FALSE, check.names = FALSE,
               colClasses = "character")
need_map <- c("feature_id", "mapped_entrez", "status")
if (!all(need_map %in% names(mp))) stop("Mapping file lacks required columns")

# Match the EXACT annotated UniProt entry used in the previous coverage audit.
# Do not combine or select among multiple proteoforms without verification.
raw_accessions <- unique(mp$feature_id[mp$mapped_entrez == "201294" &
                                      mp$status == "UNIQUE"])
raw_accessions <- intersect(raw_accessions, rownames(RAW))
if (length(raw_accessions) != 1L)
  stop("UNC13D has ", length(raw_accessions), " uniquely mapped raw accessions (expected 1). ",
       "Inspect 15_raw_uniprot_entrez_mapping.csv; do NOT sum intensities arbitrarily.")
if (!"201294" %in% rownames(PRO)) stop("Processed UNC13D Entrez ID 201294 not found")
cat("UNC13D raw UniProt accession: ", raw_accessions, "\n", sep = "")

D <- read.csv(file.path(ROOT, "donor.csv"), stringsAsFactors = FALSE,
              check.names = FALSE)
COV <- c("donorage", "donorsex", "bodymassindex",
         "predistributionculturetime", "percentieqrecoverypostculture")
required <- c("record_id", "diagnosis_computed", COV)
if (!all(required %in% names(D)))
  stop("Missing metadata: ", paste(setdiff(required, names(D)), collapse = ", "))
if (anyDuplicated(D$record_id)) stop("Duplicate donor record IDs")
D$grp <- unname(c("None" = "ND", "Type2" = "T2D")[
  trimws(as.character(D$diagnosis_computed))])
shared <- intersect(colnames(RAW), colnames(PRO))
keep <- !is.na(D$grp) & D$record_id %in% shared &
        complete.cases(D[, COV, drop = FALSE])
d <- D[keep, , drop = FALSE]
d$grp <- factor(d$grp, levels = c("ND", "T2D"))
d$donorsex <- factor(d$donorsex)
ids <- as.character(d$record_id)
if (anyNA(ids) || anyDuplicated(ids)) stop("Invalid matched donor IDs")
if (!all(c("ND", "T2D") %in% d$grp)) stop("Missing comparison group")
cat("Eligible discovery donors: ")
print(table(d$grp))
if (sum(d$grp == "ND") != 87L || sum(d$grp == "T2D") != 16L)
  stop("Expected 87 ND and 16 T2D; check metadata before continuing")

raw_y <- as.numeric(RAW[raw_accessions, ids])
pro_y <- as.numeric(PRO["201294", ids])
observed <- is.finite(raw_y) & raw_y > 0
if (any(!is.finite(pro_y)))
  warning("Processed UNC13D has nonfinite values; those are excluded from model fits")
cat("Raw UNC13D measurement availability by group:\n")
print(table(group = d$grp, raw_observed = observed))
stopifnot(sum(observed & d$grp == "ND") == 59L,
          sum(observed & d$grp == "T2D") == 15L)

# The raw matrix is on an intensity scale. Normalize as a sensitivity check:
#  (i) log2(raw UNC13D) without sample normalization
#  (ii) log2(raw UNC13D) - each donor's across-protein median + grand median
# These are NOT presumed to equal the source project's processing pipeline.
R_all <- RAW[, ids, drop = FALSE]
R_all[!is.finite(R_all) | R_all <= 0] <- NA_real_
L_all <- log2(R_all)
sample_median <- apply(L_all, 2, median, na.rm = TRUE)
sample_coverage <- colSums(is.finite(L_all))
if (any(!is.finite(sample_median))) stop("Unable to compute median for a donor")
grand_median <- median(sample_median)
raw_log <- rep(NA_real_, length(ids))
raw_log[observed] <- log2(raw_y[observed])
raw_median_norm <- raw_log - sample_median[ids] + grand_median

# Raw measurement availability: Fisher exact is UNADJUSTED and exploratory.
cont <- table(group = d$grp, raw_observed = factor(observed, levels = c(TRUE, FALSE)))
ft <- fisher.test(cont)
cat(sprintf("\nRaw missingness Fisher exact P = %.5f, OR (ND vs T2D observed odds) = %.4f\n",
            ft$p.value, unname(ft$estimate)))

# Detect candidate batch identifiers, but do not guess which is a true
# mass-spectrometry/acquisition batch. User must set BATCH_COL if documented.
batch_candidates <- grep("batch|plate|run|plex|instrument|acquis|mass.spec|ms_batch",
                         names(D), ignore.case = TRUE, value = TRUE)
cat("Potential discovery batch fields: ",
    if (length(batch_candidates)) paste(batch_candidates, collapse = ", ") else "NONE", "\n", sep = "")
batch_summary <- NULL
if (nzchar(BATCH_COL)) {
  if (!BATCH_COL %in% names(d)) stop("BATCH_COL is not a donor metadata column")
  b <- as.character(d[[BATCH_COL]])
  if (anyNA(b)) stop("Missing values in chosen batch field; review first")
  tmp <- data.frame(batch = b, group = d$grp, observed = observed)
  by_counts <- aggregate(rep(1L, nrow(tmp)),
                         tmp[, c("batch", "group", "observed"), drop = FALSE], sum)
  names(by_counts)[4] <- "n_donors"
  batch_summary <- by_counts
  write.csv(batch_summary, file.path(OUT, "22_unc13d_missingness_by_batch.csv"),
            row.names = FALSE)
  cat("Raw missingness by VERIFIED batch:\n")
  print(xtabs(~ batch + group + observed, data = tmp))
} else {
  writeLines(c("A batch-stratified test was NOT run because no verified acquisition",
               "batch column was supplied. Inspect the candidate columns printed",
               "above and set BATCH_COL only after confirming its meaning.",
               paste("Candidate columns:", paste(batch_candidates, collapse = ", "))),
             file.path(OUT, "22_unc13d_batch_not_verified.txt"))
}

# Same nuisance covariates as the discovery analyses. We compute a residual
# SD ratio and Brown-Forsythe statistic on those residuals. The latter is
# DESCRIPTIVE: standard BF P values are not a fully validated adjusted test.
model_formula <- function(group_mean = FALSE) {
  as.formula(paste("y ~", if (group_mean) "grp +" else "",
                   "donorage + donorsex + bodymassindex +",
                   "predistributionculturetime + percentieqrecoverypostculture"))
}

estimate <- function(d0, yy, group_mean = FALSE) {
  ok <- is.finite(yy) & !is.na(d0$grp)
  dd <- d0[ok, , drop = FALSE]
  dd$y <- yy[ok]
  dd$grp <- droplevels(factor(dd$grp, levels = c("ND", "T2D")))
  if (any(table(dd$grp) < 4L) || length(levels(dd$grp)) != 2L)
    stop("Too few observations in a group")
  # Keep sex factor coding explicit and check for design aliasing.
  dd$donorsex <- factor(dd$donorsex)
  fm <- model_formula(group_mean)
  XX <- model.matrix(fm, data = dd)
  if (qr(XX)$rank != ncol(XX)) stop("Model rank deficient after donor selection")
  fit <- lm(fm, data = dd)
  rr <- residuals(fit)
  nd <- dd$grp == "ND"; td <- dd$grp == "T2D"
  s_nd <- sd(rr[nd]); s_td <- sd(rr[td])
  if (!is.finite(s_nd) || s_nd <= 0) stop("Invalid ND residual SD")
  sd_ratio <- s_td / s_nd
  zz <- numeric(length(rr))
  zz[nd] <- abs(rr[nd] - median(rr[nd]))
  zz[td] <- abs(rr[td] - median(rr[td]))
  bf_p <- anova(lm(z ~ grp, data = data.frame(z = zz, grp = dd$grp)))[["Pr(>F)"]][1]
  c(n_ND = sum(nd), n_T2D = sum(td), sd_ND = s_nd, sd_T2D = s_td,
    sd_ratio = sd_ratio, brown_forsythe_approx_p = bf_p)
}

# Cases: full processed; processed on raw-observed donors; raw log2 intensity
# and normalized raw log2 intensity on exactly the same observed donors.
cases <- list(
  list(name = "processed_all_103", data = d, y = pro_y, group_mean = FALSE),
  list(name = "processed_raw_observed_74", data = d[observed, , drop = FALSE],
       y = pro_y[observed], group_mean = FALSE),
  list(name = "raw_log2_observed_74", data = d[observed, , drop = FALSE],
       y = raw_log[observed], group_mean = FALSE),
  list(name = "raw_median_normalized_74", data = d[observed, , drop = FALSE],
       y = raw_median_norm[observed], group_mean = FALSE),
  list(name = "processed_raw_observed_74_group_mean", data = d[observed, , drop = FALSE],
       y = pro_y[observed], group_mean = TRUE),
  list(name = "raw_normalized_74_group_mean", data = d[observed, , drop = FALSE],
       y = raw_median_norm[observed], group_mean = TRUE)
)

boot_ratio <- function(d0, yy, group_mean, B = NBOOT) {
  nd <- which(d0$grp == "ND"); td <- which(d0$grp == "T2D")
  out <- rep(NA_real_, B)
  for (i in seq_len(B)) {
    ix <- c(sample(nd, length(nd), replace = TRUE),
            sample(td, length(td), replace = TRUE))
    v <- try(estimate(d0[ix, , drop = FALSE], yy[ix], group_mean)[["sd_ratio"]],
             silent = TRUE)
    if (!inherits(v, "try-error") && is.finite(v)) out[i] <- v
  }
  good <- out[is.finite(out)]
  if (length(good) < 0.9 * B) stop("Over 10% of bootstrap fits failed; inspect model")
  c(ci_low = unname(quantile(good, 0.025)),
    ci_high = unname(quantile(good, 0.975)),
    n_boot_valid = length(good), n_boot_failed = B - length(good))
}

results <- do.call(rbind, lapply(cases, function(z) {
  est <- estimate(z$data, z$y, z$group_mean)
  bt <- boot_ratio(z$data, z$y, z$group_mean)
  data.frame(analysis = z$name, group_in_location_model = z$group_mean,
             t(c(est, bt)), row.names = NULL, check.names = FALSE)
}))
for (nm in setdiff(names(results), c("analysis", "group_in_location_model"))) {
  results[[nm]] <- as.numeric(results[[nm]])
}
cat("\n===== UNC13D SENSITIVITY RESULTS =====\n")
print(results, row.names = FALSE, digits = 4)
write.csv(results, file.path(OUT, "21_unc13d_variance_sensitivity.csv"), row.names = FALSE)

# Isolate donor selection while holding the model fit fixed: compute residuals
# on all 103 processed samples, then restrict to the 74 raw-observed donors.
full_model <- lm(model_formula(FALSE), data = transform(d, y = pro_y))
rr <- residuals(full_model)
full_fit_restricted <- sd(rr[observed & d$grp == "T2D"]) /
                       sd(rr[observed & d$grp == "ND"])
cat(sprintf(paste0("\nProcessed residual SD ratio (full 103-donor nuisance fit,\n",
                   "  restricted to the 74 raw-observed donors): %.4f\n"),
            full_fit_restricted))
write.csv(data.frame(analysis = "processed_full_fit_residuals_restricted_to_74",
                     sd_ratio = full_fit_restricted),
          file.path(OUT, "23_unc13d_fixed_fit_restriction.csv"), row.names = FALSE)

# QC file: all eligible donors, including those whose raw UNC13D is missing.
qc <- d[, c("record_id", "grp", COV), drop = FALSE]
qc$raw_unc13d_observed <- observed
qc$raw_unc13d_intensity <- raw_y
qc$raw_unc13d_log2 <- raw_log
qc$raw_median_normalized_log2 <- raw_median_norm
qc$processed_unc13d <- pro_y
qc$raw_global_log2_median <- unname(sample_median[ids])
qc$raw_number_features_observed <- unname(sample_coverage[ids])
if (nzchar(BATCH_COL)) qc$verified_batch <- d[[BATCH_COL]]
write.csv(qc, file.path(OUT, "20_unc13d_donor_measurement_audit.csv"), row.names = FALSE)

writeLines(c(
  "Analysis completed. Interpretation cautions:",
  "1. The unprocessed file gives intensity, not proof of protein presence/absence.",
  "2. Raw median normalization is a sensitivity analysis, not the verified source pipeline.",
  "3. Brown-Forsythe P values on regression residuals are approximate; not confirmatory.",
  "4. Bootstrap percentile intervals refit the nuisance model but may undercover",
  "   with small T2D samples or informative missingness; use as exploratory.",
  "5. Without a verified proteomics acquisition batch, batch missingness",
  "   cannot be directly evaluated. Do not substitute RNA library batch.",
  "6. Complete-case analyses can be biased if missingness is informative.",
  "7. If the adjusted raw and processed ratios disagree, first verify source",
  "   normalization, imputation and identifier handling before interpretation."
), file.path(OUT, "24_unc13d_analysis_cautions.txt"))
cat("\nWrote files 20, 21, 23 and 24 in: ", OUT, "\n", sep = "")
