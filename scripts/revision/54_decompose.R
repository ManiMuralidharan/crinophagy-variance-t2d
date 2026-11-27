# =====================================================================
# 54_decompose.R
#
# Why does UNC13D have an SD ratio of 2.00 but a dispersion-model
# t of -0.12? Three candidate explanations, separated here:
#
#   (A) composition genuinely explains it
#   (B) over-adjustment: group is near-aliased with imm/beta
#   (C) the SD ratio is driven by a few extreme donors, and a robust
#       statistic does not see it
#
# Diagnostics:
#   1. how collinear group is with imm and beta, within batch
#   2. nested dispersion models, so you can see where t collapses
#   3. SD ratio vs robust MAD ratio for each target gene
#   4. the actual residual values for UNC13D, donor by donor
#
# RUN:  source("40_common.R"); source("54_decompose.R")
# =====================================================================

stopifnot(exists("mat"), exists("pd2"), exists("crino"))
set.seed(1)

grp  <- factor(as.character(pd2$group), levels = c("ND", "T2D"))
SAMP <- rownames(pd2)

TARG <- c(UNC13D = "ENSG00000092929", STX7 = "ENSG00000079950",
          GABARAP = "ENSG00000170296", CTSD = "ENSG00000117984")
TARG <- TARG[TARG %in% rownames(mat)]
cat("targets resolved:", paste(names(TARG), collapse = ", "), "\n")

col_med <- function(M) apply(M, 2, median)
spread_matrix <- function(RES, lab) {
  Z <- RES
  for (g in unique(lab)) {
    ix <- lab == g
    Z[ix, ] <- abs(sweep(RES[ix, , drop = FALSE], 2,
                         col_med(RES[ix, , drop = FALSE])))
  }
  Z
}
coef_t <- function(X, Y, j) {
  f <- lm.fit(X, Y); n <- nrow(X); p <- ncol(X)
  s2 <- colSums(f$residuals^2) / (n - p)
  XtXi <- chol2inv(chol(crossprod(X)))
  f$coefficients[j, ] / sqrt(s2 * XtXi[j, j])
}


## ---- 1. How aliased is group with composition? ----------------------

cat("\n===== 1. COLLINEARITY: group vs composition =====\n")
m_imm  <- lm(imm  ~ grp + batch, data = pd2)
m_beta <- lm(beta ~ grp + batch, data = pd2)
cat(sprintf("  imm  ~ grp + batch : R2 = %.3f, grp coefficient P = %.4g\n",
            summary(m_imm)$r.squared,  summary(m_imm)$coefficients["grpT2D", 4]))
cat(sprintf("  beta ~ grp + batch : R2 = %.3f, grp coefficient P = %.4g\n",
            summary(m_beta)$r.squared, summary(m_beta)$coefficients["grpT2D", 4]))

Xfull <- model.matrix(~ grp + imm + beta + batch, data = pd2)
r2_grp <- summary(lm(Xfull[, "grpT2D"] ~ Xfull[, setdiff(colnames(Xfull),
                      c("(Intercept)", "grpT2D"))]))$r.squared
cat(sprintf("\n  R2 of grpT2D regressed on all other model terms: %.3f\n", r2_grp))
cat(sprintf("  variance inflation factor for group: %.2f\n", 1 / (1 - r2_grp)))
cat("  VIF above ~5 means the group term is substantially aliased and the\n")
cat("  fully adjusted dispersion model is over-adjusted. Above ~10, severely.\n")


## ---- 2. Nested dispersion models ------------------------------------

RESL <- residual_matrix(~ batch + imm + beta)
Z    <- spread_matrix(RESL, as.character(grp))

models <- list(
  "z ~ group"                      = ~ grp,
  "z ~ group + batch"              = ~ grp + batch,
  "z ~ group + batch + imm"        = ~ grp + batch + imm,
  "z ~ group + batch + beta"       = ~ grp + batch + beta,
  "z ~ group + batch + imm + beta" = ~ grp + batch + imm + beta
)

cat("\n===== 2. NESTED DISPERSION MODELS: t on the group term =====\n")
res_tab <- lapply(names(models), function(nm) {
  X <- model.matrix(models[[nm]], data = pd2)
  j <- which(colnames(X) == "grpT2D")
  tv <- coef_t(X, Z, j)
  c(setNames(sapply(TARG, function(e) tv[e]), names(TARG)),
    crino_median = median(tv[match(crino, names(tv))], na.rm = TRUE),
    genomewide_median = median(tv, na.rm = TRUE),
    UNC13D_pct = 100 * mean(tv <= tv[TARG[["UNC13D"]]], na.rm = TRUE))
})
res_tab <- do.call(rbind, res_tab)
rownames(res_tab) <- names(models)
print(round(res_tab, 3))
cat("\n  Read down the UNC13D column. If t is large in the first rows and\n")
cat("  collapses only when imm/beta enter, that is explanation (A) or (B),\n")
cat("  and the VIF above tells you which. If UNC13D is small even in\n")
cat("  'z ~ group', composition is not the cause and explanation (C) is\n")
cat("  in play.\n")
cat("  genomewide_median is the calibration check: it should sit near 0.\n")


## ---- 3. SD ratio vs robust ratio ------------------------------------

cat("\n===== 3. SD ratio vs robust (MAD) ratio =====\n")
rob <- do.call(rbind, lapply(names(TARG), function(s) {
  e <- TARG[[s]]
  r <- RESL[, e]
  a <- r[grp == "T2D"]; b <- r[grp == "ND"]
  # trimmed SD: drop the two most extreme residuals in each group
  tsd <- function(v) { v <- sort(abs(v - median(v))); sd(v[seq_len(length(v) - 2)]) }
  data.frame(gene = s,
             sd_ratio    = sd(a) / sd(b),
             mad_ratio   = mad(a) / mad(b),
             iqr_ratio   = IQR(a) / IQR(b),
             trimmed_ratio = tsd(a) / tsd(b),
             n_T2D = length(a), n_ND = length(b))
}))
print(rob, row.names = FALSE, digits = 4)
cat("\n  A large sd_ratio with mad_ratio and iqr_ratio near 1 means the\n")
cat("  signal lives in a few extreme donors, not in a broad difference in\n")
cat("  spread. That is explanation (C), and it is the one Referee 2's\n")
cat("  minor comment 4 was circling.\n")


## ---- 4. UNC13D donor by donor ---------------------------------------

cat("\n===== 4. UNC13D residuals, donor by donor =====\n")
r <- RESL[, TARG[["UNC13D"]]]
dd <- data.frame(donor = SAMP, group = as.character(grp), residual = r,
                 imm = pd2$imm, beta = pd2$beta, batch = as.character(pd2$batch))
dd$abs_dev <- abs(dd$residual - ave(dd$residual, dd$group, FUN = median))
dd <- dd[order(-dd$abs_dev), ]
cat("  ten largest absolute deviations:\n")
print(head(dd, 10), row.names = FALSE, digits = 3)

a <- r[grp == "T2D"]; b <- r[grp == "ND"]
cat(sprintf("\n  T2D: n=%d sd=%.3f mad=%.3f   ND: n=%d sd=%.3f mad=%.3f\n",
            length(a), sd(a), mad(a), length(b), sd(b), mad(b)))
drop2 <- function(v) { o <- order(abs(v - median(v)), decreasing = TRUE); v[-o[1:2]] }
cat(sprintf("  after removing the 2 most extreme T2D donors: sd ratio %.3f (was %.3f)\n",
            sd(drop2(a)) / sd(b), sd(a) / sd(b)))

write.csv(as.data.frame(res_tab), "C_nested_dispersion.csv")
write.csv(rob, "C_robust_ratios.csv", row.names = FALSE)
write.csv(dd,  "C_unc13d_donor_level.csv", row.names = FALSE)

cat("\n===== WHAT EACH PATTERN MEANS FOR THE PAPER =====\n")
cat("  (A) composition explains it: UNC13D t large without imm/beta, small\n")
cat("      with them, and VIF low. The referees' concern is correct and the\n")
cat("      gene-level claim needs heavy qualification.\n")
cat("  (B) over-adjustment: same collapse but VIF high. Then the fully\n")
cat("      adjusted dispersion model is not the right test, and the\n")
cat("      published Brown-Forsythe result stands with a stated caveat.\n")
cat("  (C) outlier-driven: sd_ratio >> mad_ratio, and dropping two donors\n")
cat("      moves the ratio a lot. Then report a robust estimate alongside\n")
cat("      the SD ratio and say plainly how many donors carry the signal.\n")
cat("  These are not mutually exclusive. Send me the whole output.\n")
