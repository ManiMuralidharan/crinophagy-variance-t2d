# 72_build_integrated_discovery_supplement.R
# IJMS revision: reproducible, full-precision merger of HISTORICAL processed
# 26-protein testing and NEW exploratory raw 31-gene sensitivity analysis.
# This script reads original CSV files and never overwrites them.
# Resulting P values on regression residuals are not a fully calibrated
# covariate-adjusted test; do NOT combine P values across methods.

ROOT <- "."
OUT <- file.path(ROOT, "reviewer_validation_output")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
HIST_FILE <- file.path(ROOT, "per_protein_variance.csv")
RAW_FILE <- file.path(OUT, "41_raw31_variance_and_BH_sensitivity.csv")
VAR_FILE <- file.path(OUT, "42_raw31_feature_variance_rank.csv")
COV_FILE <- file.path(OUT, "16_raw_and_processed_panel_coverage.csv")
files <- c(HIST_FILE, RAW_FILE, VAR_FILE, COV_FILE)
if (any(!file.exists(files))) stop("Missing required file(s):\n",paste(files[!file.exists(files)],collapse="\n"))

H <- read.csv(HIST_FILE, stringsAsFactors=FALSE, check.names=FALSE)
R <- read.csv(RAW_FILE, stringsAsFactors=FALSE, check.names=FALSE)
V <- read.csv(VAR_FILE, stringsAsFactors=FALSE, check.names=FALSE)
C <- read.csv(COV_FILE, stringsAsFactors=FALSE, check.names=FALSE)
stopifnot(all(c("symbol","set","sd_ratio","p","q") %in% names(H)))
stopifnot(all(c("gene","entrez","raw_accessions","n_raw_accessions","processed_rows","n_model","n_ND","n_preT2D","n_T2D","sd_ratio","bf_p_approx","q_BH_evaluable","sd_ratio_group_mean") %in% names(R)))
stopifnot(all(c("gene","accession","n_positive_finite_all_raw_donors","variance_percentile_from_low","at_or_below_20pct_reference") %in% names(V)))
stopifnot(all(c("gene","set","raw_nonmissing_ND","raw_nonmissing_T2D","processed_rows") %in% names(C)))
H <- H[H$set == "crinophagy", , drop=FALSE]
C <- C[C$set == "crinophagy", , drop=FALSE]
if (nrow(H)!=26L || nrow(R)!=31L || nrow(C)!=31L)
 stop(sprintf("Unexpected gene counts: historical=%d, raw=%d, coverage=%d",nrow(H),nrow(R),nrow(C)))
if (anyDuplicated(H$symbol) || anyDuplicated(R$gene) || anyDuplicated(C$gene))
 stop("Duplicate gene names in source tables")
if (!identical(sort(R$gene),sort(C$gene))) stop("Raw and frozen coverage gene sets do not match")
if (!setequal(H$symbol, R$gene[R$processed_rows > 0]))
 stop("Historical tested-gene list does not match processed features; audit before publication")
if (anyNA(R$bf_p_approx) || anyNA(R$q_BH_evaluable))
 stop("Expected 31 evaluable raw tests; inspect untestable genes before proceeding")

# Verify the original historical FDR family and the NEW raw 31-gene family.
check_bh <- function(p,q,label) {
  if (any(!is.finite(p)) || any(!is.finite(q))) stop("Invalid P/q values in ",label)
  expected <- p.adjust(p, method="BH")
  delta <- max(abs(expected-q))
  cat(sprintf("%s: n=%d, max(abs(BH(recomputed)-supplied))=%.8g\n",label,length(p),delta))
  if (delta > 2e-6) stop("BH discrepancy in ",label,"; do not merge until resolved")
}
check_bh(H$p,H$q,"Historical processed crinophagy")
check_bh(R$bf_p_approx,R$q_BH_evaluable,"Exploratory raw crinophagy")

match_h <- match(R$gene,H$symbol)
match_c <- match(R$gene,C$gene)
merged <- data.frame(
 gene=R$gene,
 historical_processed_tested=!is.na(match_h),
 historical_fit_n=ifelse(!is.na(match_h),123L,NA_integer_),
 historical_SD_ND_n=ifelse(!is.na(match_h),87L,NA_integer_),
 historical_SD_T2D_n=ifelse(!is.na(match_h),16L,NA_integer_),
 historical_processed_SD_ratio=H$sd_ratio[match_h],
 historical_processed_BF_P=H$p[match_h],
 historical_processed_BH_q_26=H$q[match_h],
 raw_n_model=R$n_model,
 raw_ND_n=R$n_ND,
 raw_preT2D_n=R$n_preT2D,
 raw_T2D_n=R$n_T2D,
 raw_log2_median_normalized_SD_ratio=R$sd_ratio,
 raw_BF_approx_P=R$bf_p_approx,
 raw_BH_approx_q_31=R$q_BH_evaluable,
 raw_group_mean_SD_ratio=R$sd_ratio_group_mean,
 raw_mapped_accessions=R$raw_accessions,
 raw_accession_count=R$n_raw_accessions,
 raw_nonmissing_ND_n=C$raw_nonmissing_ND[match_c],
 raw_nonmissing_T2D_n=C$raw_nonmissing_T2D[match_c],
 processed_rows=R$processed_rows,
 stringsAsFactors=FALSE
)
merged$note <- ifelse(merged$historical_processed_tested,
                     "Historical processed-matrix test; raw values are separate exploratory sensitivity",
                     "Not tested historically; feature absent from processed matrix")
write.csv(merged,file.path(OUT,"45_combined_processed26_raw31_discovery.csv"),row.names=FALSE,na="")

# Preserve accessions and variance-rank details of the five omitted genes.
missing <- merged$gene[!merged$historical_processed_tested]
if (length(missing)!=5L) stop("Expected exactly five historical exclusions")
vv <- V[V$gene %in% missing,,drop=FALSE]
if (!all(missing %in% vv$gene)) stop("Some excluded genes lack raw variance-rank records")
vv <- vv[match(missing,vv$gene),,drop=FALSE]
excl <- data.frame(
 gene=vv$gene,
 accession=vv$accession,
 raw_positive_finite_all_donors=vv$n_positive_finite_all_raw_donors,
 raw_variance_percentile_from_low=vv$variance_percentile_from_low,
 at_or_below_bottom20_raw=vv$at_or_below_20pct_reference,
 processed_rows=merged$processed_rows[match(vv$gene,merged$gene)],
 inference="Consistent with low-variance filtering; exact provider exclusion mechanism not confirmed",
 stringsAsFactors=FALSE
)
write.csv(excl,file.path(OUT,"46_five_excluded_processing_audit.csv"),row.names=FALSE)

writeLines(c(
 "REVIEWER 2 METHODS/RESULTS NOTES (Script 72)",
 "Historical processed: proc_prot.rds, 26 evaluated genes, nuisance fit n=123; SD ratio uses 87 ND/16 T2D.",
 "Raw sensitivity: unproc_prot.rds, log2 + donor-median centered, no imputation, 31/31 genes evaluable.",
 "Raw nuisance fit uses all available ND/preT2D/T2D donors with complete clinical covariates; SD ratio compares ND/T2D only.",
 "Historical original Brown-Forsythe P values and q values are copied from per_protein_variance.csv and BH-verified.",
 "Raw Brown-Forsythe residual P values are APPROXIMATE; q values correct these approximate P values within the 31-gene raw panel.",
 "NO inference that a raw q below 0.05 is a calibrated confirmatory test.",
 "The five processed-excluded genes have raw-positive signals in all 129 donors, but low raw variance rank is only CONSISTENT with provider filtering.",
 "Do not label genes with absent processed rows as undetected or negative.",
 "Patient-level batch, original export/source version, residual-leverage testing and bootstrap coverage remain unresolved.",
 "UNC13D processed q 0.02154 versus raw 31-panel q 0.11855; do not describe it as robustly FDR significant across pipelines.",
 "Original historical source files were not modified."
),file.path(OUT,"47_integrated_discovery_analysis_notes.txt"))
cat("\n===== KEY COMPARISON =====\n")
print(merged[merged$gene %in% c("UNC13D","RAB3A"),
 c("gene","historical_processed_SD_ratio","historical_processed_BF_P","historical_processed_BH_q_26",
   "raw_ND_n","raw_T2D_n","raw_log2_median_normalized_SD_ratio","raw_BF_approx_P","raw_BH_approx_q_31")],row.names=FALSE,digits=7)
cat("\n===== FIVE EXCLUSIONS =====\n")
print(excl,row.names=FALSE)
cat("\nCreated CSVs 45, 46 and notes 47 in: ",OUT,"\n",sep="")
