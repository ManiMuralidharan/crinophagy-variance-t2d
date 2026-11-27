# =====================================================================
# 53_dispersion_model.R
#
# Tests disease association in DISPERSION directly, with composition
# adjusted in the dispersion model rather than only in the location
# model. This avoids the exchangeability question entirely: no labels
# are permuted, so no null of exchangeable labels is required.
#
# For each gene:
#   1. location model   y ~ batch + imm + beta      -> residuals r
#   2. spread variable  z = |r - median(r) within group|   (Brown-Forsythe)
#   3. dispersion model z ~ group + imm + beta + batch
#      The coefficient on `group` is the covariate-adjusted
#      disease-dispersion association. Its t statistic is the gene's
#      evidence.
#
# Set-level questions then need no permutation:
#   - are crinophagy t statistics shifted relative to all other genes?
#   - are they shifted relative to the comparator set?
#   Both by Wilcoxon rank-sum, which conditions on nothing.
#
# Also runs: group-aware location model (recommendation 1), and
# leave-one-donor-out influence (recommendation 3).
#
# RUN:  source("40_common.R"); source("53_dispersion_model.R")
# =====================================================================

stopifnot(exists("mat"), exists("pd2"), exists("crino"), exists("macro"))
set.seed(1)

SAMP <- rownames(pd2)
grp  <- factor(as.character(pd2$group), levels = c("ND", "T2D"))

TARGETS <- c(UNC13D = "ENSG00000092929", STX7 = "ENSG00000079950")
GAB     <- "GABARAP"   # resolved below if possible


## ---- helpers --------------------------------------------------------

colsd <- function(M) sqrt(colSums(sweep(M, 2, colMeans(M))^2) / (nrow(M) - 1))

col_med <- function(M) apply(M, 2, median)

spread_matrix <- function(RES, lab) {
  # |residual - group median|, computed per gene within each group
  Z <- RES
  for (g in unique(lab)) {
    ix <- lab == g
    Z[ix, ] <- abs(sweep(RES[ix, , drop = FALSE], 2, col_med(RES[ix, , drop = FALSE])))
  }
  Z
}

# vectorised t statistic for one named coefficient across all genes
coef_t <- function(X, Y, which_col) {
  fit  <- lm.fit(X, Y)
  n    <- nrow(X); p <- ncol(X)
  rss  <- colSums(fit$residuals^2)
  s2   <- rss / (n - p)
  XtXi <- chol2inv(chol(crossprod(X)))
  se   <- sqrt(s2 * XtXi[which_col, which_col])
  list(est = fit$coefficients[which_col, ], se = se,
       t = fit$coefficients[which_col, ] / se, df = n - p)
}


## ---- 1. Dispersion model, composition adjusted in BOTH stages -------

RESL <- residual_matrix(~ batch + imm + beta)      # location (from 40_common.R)
Z    <- spread_matrix(RESL, as.character(grp))

Xd   <- model.matrix(~ grp + imm + beta + batch, data = pd2)
gcol <- which(colnames(Xd) == "grpT2D")
stopifnot(length(gcol) == 1)

fit  <- coef_t(Xd, Z, gcol)
tval <- fit$t
pval <- 2 * pt(abs(tval), df = fit$df, lower.tail = FALSE)
cat(sprintf("Dispersion model fitted for %d genes, df = %d\n", length(tval), fit$df))


## ---- 2. Set-level, no permutation required --------------------------

ci <- match(crino, names(tval)); ci <- ci[!is.na(ci)]
mi <- match(macro, names(tval)); mi <- mi[!is.na(mi)]
rest <- setdiff(seq_along(tval), c(ci, mi))

cat("\n===== SET-LEVEL: covariate-adjusted dispersion association =====\n")
cat(sprintf("  crinophagy (n=%d)  median t = %+.3f\n", length(ci), median(tval[ci])))
cat(sprintf("  comparator (n=%d)  median t = %+.3f\n", length(mi), median(tval[mi])))
cat(sprintf("  all other genes    median t = %+.3f\n", median(tval[rest])))

w1 <- wilcox.test(tval[ci], tval[rest], alternative = "greater")
w2 <- wilcox.test(tval[ci], tval[mi],   alternative = "greater")
w3 <- wilcox.test(tval[mi], tval[rest], alternative = "greater")
cat(sprintf("\n  crinophagy vs all other genes : Wilcoxon one-sided P = %.4g\n", w1$p.value))
cat(sprintf("  crinophagy vs comparator      : Wilcoxon one-sided P = %.4g\n", w2$p.value))
cat(sprintf("  comparator vs all other genes : Wilcoxon one-sided P = %.4g\n", w3$p.value))
cat("\n  These test whether the SET's disease-dispersion statistics are\n")
cat("  shifted. No labels are permuted, so exchangeability is not assumed.\n")


## ---- 3. Gene level ---------------------------------------------------

want <- c(TARGETS, setNames(
  names(tval)[match(GAB, names(tval))], GAB))     # GABARAP by Ensembl if keyed so
want <- want[!is.na(want) & want %in% names(tval)]

# add GABARAP by its Ensembl id if the symbol lookup failed
GAB_ENS <- "ENSG00000170296"
if (!GAB %in% names(want) && GAB_ENS %in% names(tval))
  want <- c(want, setNames(GAB_ENS, "GABARAP"))

cat("\n===== GENE LEVEL (covariate-adjusted dispersion) =====\n")
gl <- do.call(rbind, lapply(names(want), function(s) {
  e <- want[[s]]
  data.frame(gene = s, ensembl = e,
             est = unname(fit$est[e]), se = unname(fit$se[e]),
             t = unname(tval[e]), p = unname(pval[e]),
             pct_t = 100 * mean(tval <= tval[e]))
}))
print(gl, row.names = FALSE, digits = 4)
cat("\n  pct_t is the gene's percentile among all genes on this statistic.\n")


## ---- 4. Recommendation 1: group-aware location model ----------------
# If group is included in the location model, group mean differences
# cannot leak into the spread variable at all.

RESL_g <- residual_matrix(~ grp + batch + imm + beta)
Zg     <- spread_matrix(RESL_g, as.character(grp))
fit_g  <- coef_t(Xd, Zg, gcol)

cat("\n===== SENSITIVITY: group included in the location model =====\n")
cat(sprintf("  crinophagy median t : %+.3f  (was %+.3f)\n",
            median(fit_g$t[ci]), median(tval[ci])))
cat(sprintf("  comparator median t : %+.3f  (was %+.3f)\n",
            median(fit_g$t[mi]), median(tval[mi])))
cat(sprintf("  crinophagy vs rest  : P = %.4g  (was %.4g)\n",
            wilcox.test(fit_g$t[ci], fit_g$t[-c(ci, mi)],
                        alternative = "greater")$p.value, w1$p.value))
for (s in names(want))
  cat(sprintf("  %-8s t = %+.3f  (was %+.3f)\n", s,
              fit_g$t[want[[s]]], tval[want[[s]]]))


## ---- 5. Recommendation 3: leave-one-donor-out -----------------------

cat("\n===== INFLUENCE: leave one donor out =====\n")
loo <- do.call(rbind, lapply(seq_len(nrow(pd2)), function(i) {
  keep <- setdiff(seq_len(nrow(pd2)), i)
  pdx  <- droplevels(pd2[keep, ])
  Xl   <- model.matrix(~ batch + imm + beta, data = pdx)
  Yl   <- t(mat[, SAMP[keep], drop = FALSE])
  Rl   <- lm.fit(Xl, Yl)$residuals
  Zl   <- spread_matrix(Rl, as.character(pdx$group))
  Xdl  <- model.matrix(~ group + imm + beta + batch, data = pdx)
  gc2  <- which(colnames(Xdl) == "groupT2D")
  fl   <- coef_t(Xdl, Zl, gc2)
  cil  <- match(crino, names(fl$t)); cil <- cil[!is.na(cil)]
  data.frame(dropped = SAMP[i], group = as.character(pd2$group)[i],
             crino_median_t = median(fl$t[cil]),
             UNC13D_t = unname(fl$t[TARGETS[["UNC13D"]]]))
}))
loo <- loo[order(loo$UNC13D_t), ]
cat("  five donors whose removal lowers UNC13D most:\n")
print(head(loo, 5), row.names = FALSE, digits = 4)
cat("  five whose removal raises it most:\n")
print(tail(loo, 5), row.names = FALSE, digits = 4)
cat(sprintf("\n  full-data UNC13D t = %+.3f; range across leave-one-out: %+.3f to %+.3f\n",
            tval[TARGETS[["UNC13D"]]], min(loo$UNC13D_t), max(loo$UNC13D_t)))
cat(sprintf("  crinophagy median t range across leave-one-out: %+.3f to %+.3f\n",
            min(loo$crino_median_t), max(loo$crino_median_t)))

save(tval, pval, fit, fit_g, loo, file = "C_dispersion_model.RData")
write.csv(gl,  "C_dispersion_gene_level.csv", row.names = FALSE)
write.csv(loo, "C_dispersion_loo.csv", row.names = FALSE)

cat("\n===== HOW TO READ THIS =====\n")
cat("  This is the test the set-level claim should stand or fall on,\n")
cat("  because composition is adjusted where it acts (the dispersion) and\n")
cat("  no exchangeability assumption is made.\n")
cat("  If crinophagy vs all-other-genes is significant, the set-level claim\n")
cat("  survives in a defensible form and the permutation result becomes a\n")
cat("  reported sensitivity analysis.\n")
cat("  If it is not, report both and reposition the set-level finding as\n")
cat("  exploratory, with UNC13D carrying the paper.\n")
cat("\nWrote C_dispersion_gene_level.csv, C_dispersion_loo.csv, C_dispersion_model.RData\n")
