# 45_validate_C1_C2.R
# Independent diagnostic for label-permutation set statistics.
# Place beside 40_common.R and session_full.RData, then run:
# Rscript 45_validate_C1_C2.R
# It does NOT change the original data or analysis files.

source("40_common.R")

cat("\n==== 1. Batch-by-status table ====\n")
tab <- table(batch = pd2$batch,
             status = factor(as.character(pd2$group), levels = c("ND", "T2D")))
print(tab)
cat("Total donors:", nrow(pd2), "\n")
cat("Mixed batches:", sum(rowSums(tab > 0) == 2), "of", nrow(tab), "\n")
cat("Donors in mixed batches:", sum(rowSums(tab)[rowSums(tab > 0) == 2]), "\n")
cat("log10(number of distinct within-batch group allocations):",
    sum(lchoose(rowSums(tab), tab[, "T2D"])) / log(10), "\n")

cat("\n==== 2. Check gene-set membership ====\n")
if (!is.character(crino) || !is.character(macro)) {
  stop("Expected character gene-name vectors 'crino' and 'macro'; inspect str(crino) and str(macro).")
}
cat("Crinophagy gene names:", length(crino), "unique:", length(unique(crino)), "\n")
cat("Comparator gene names:", length(macro), "unique:", length(unique(macro)), "\n")
cat("Overlap between gene sets:", length(intersect(crino, macro)), "\n")
missing_genes <- setdiff(unique(c(crino, macro)), colnames(RES2))
if (length(missing_genes) > 0) stop("Panel genes not in residual matrix: ",
                                    paste(missing_genes, collapse = ", "))
panel <- unique(c(crino, macro))
stopifnot(nrow(RES2) == nrow(pd2),
          identical(as.character(pd2$group), OBS_LAB),
          !anyNA(RES2[, panel, drop = FALSE]))

cat("\n==== 3. Location-model rank and leverage ====\n")
X <- model.matrix(M2, data = pd2)
qX <- qr(X)
cat("Design rank:", qX$rank, "/", ncol(X), "columns\n")
if (qX$rank < ncol(X)) warning("Rank-deficient location model. Covariates are aliased.")
hat_diag <- rowSums(qr.Q(qX)[, seq_len(qX$rank), drop = FALSE]^2)
print(summary(hat_diag))
cat("5 highest leverage donors:\n")
print(sort(setNames(hat_diag, rownames(pd2)), decreasing = TRUE)[seq_len(min(5, length(hat_diag)))])

cat("\n==== 4. Composition-score balance ====\n")
for (v in c("imm", "beta")) {
  cat("\n", v, "by disease status:\n")
  print(by(pd2[[v]], pd2$group,
           function(x) c(n = length(x), mean = mean(x), sd = sd(x), median = median(x))))
  cat("\n", v, "within batch and status (mean):\n")
  print(with(pd2, tapply(pd2[[v]], list(pd2$batch, pd2$group), mean)))
}

cat("\n==== 5. Recompute observed set statistics ====\n")
R_panel <- RES2[, panel, drop = FALSE]
get_stat <- function(lab) {
  stopifnot(length(lab) == nrow(R_panel))
  sr <- sd_ratios(R_panel, lab)
  if (any(!is.finite(sr))) stop("Non-finite SD ratio detected")
  a <- median(sr[crino])
  b <- median(sr[macro])
  c(crinophagy = a, comparator = b, difference = a - b)
}
obs <- get_stat(OBS_LAB)
print(obs)
cat("Compare to expected: crinophagy 1.1783, comparator 1.0365, difference 0.1418.\n")

cat("\n==== 6. Independent within-batch permutation diagnostic ====\n")
# Set B to match the production script when checking exact Monte Carlo values.
B <- 9999L
set.seed(1001)
perm <- vapply(seq_len(B), function(i) get_stat(permute_labels()),
               FUN.VALUE = c(crinophagy = 0, comparator = 0, difference = 0))
null_medians <- apply(perm, 1, median)
null_ranges <- t(apply(perm, 1, quantile, probs = c(0.025, 0.975)))
p_upper <- (1 + rowSums(perm >= obs)) / (B + 1)
cat("Observed:\n"); print(obs)
cat("Permutation null medians:\n"); print(null_medians)
cat("Permutation 95% null ranges:\n"); print(null_ranges)
cat("One-sided upper-tail randomization P values:\n"); print(p_upper)
# Centred absolute-tail calculation is descriptive only. A global-sharp-null
# disease-label randomization is NOT automatically a valid test of the
# composite hypothesis that the two sets have equal NONZERO disease effects.
null_difference_center <- median(perm["difference", ])
p_difference_centered <-
  (1 + sum(abs(perm["difference", ] - null_difference_center) >=
           abs(obs["difference"] - null_difference_center))) / (B + 1)
cat("Centered two-sided difference P (global-null randomization):",
    p_difference_centered, "\n")
cat("IMPORTANT: A 95% permutation null range is not an effect-size confidence interval.\n")

cat("\n==== 7. Selected gene SD ratios (not P values) ====\n")
gene_ratios <- sd_ratios(RES2[, panel, drop = FALSE], OBS_LAB)
for (g in c("UNC13D", "STX7", "RAB3A", "GABARAP")) {
  if (g %in% names(gene_ratios)) cat(g, sprintf("%.4f", gene_ratios[g]), "\n")
}
cat("\nInterpret against the batch/composition diagnostics; do NOT use a nonsignificant P value as proof of equal variances.\n")
