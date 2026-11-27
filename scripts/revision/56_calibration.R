# =====================================================================
# 56_calibration.R      RESOLVES R2-M2 "calibration unresolved"
#
# You have two tests disagreeing on UNC13D:
#   Brown-Forsythe on residuals          P = 0.028
#   dispersion model with imm + beta     P = 0.90
# The audit correctly refuses to call either one right. This settles it
# by simulation: generate data under a KNOWN null on your actual design
# and measure each test's type-I error rate.
#
# Scenario A  no group effect, dispersion independent of composition
#             -> both tests should give ~5% at alpha 0.05
# Scenario B  no group effect, but dispersion DEPENDS on beta-cell score
#             -> this is your real situation. Whichever test still gives
#                ~5% is the trustworthy one; an inflated rate means that
#                test manufactures significance from composition alone.
# Scenario C  a real group dispersion effect, no composition dependence
#             -> power check, so an under-powered test is not mistaken
#                for a well-calibrated one.
#
# RUN:  source("40_common.R"); source("56_calibration.R")
# Time: a few minutes.
# =====================================================================

stopifnot(exists("pd2"), exists("mat"))
set.seed(1)

NGENE <- 4000        # simulated genes per scenario
ALPHA <- 0.05

grp  <- factor(as.character(pd2$group), levels = c("ND", "T2D"))
n    <- nrow(pd2)
is_t <- grp == "T2D"

# standardised composition scores, so the effect sizes below are interpretable
zb <- as.numeric(scale(pd2$beta))
zi <- as.numeric(scale(pd2$imm))

Xloc <- model.matrix(~ batch + imm + beta, data = pd2)          # location model
Xdis <- model.matrix(~ grp + imm + beta + batch, data = pd2)    # dispersion model
jdis <- which(colnames(Xdis) == "grpT2D")

col_med <- function(M) apply(M, 2, median)

spread <- function(R) {
  Z <- R
  for (g in levels(grp)) {
    ix <- grp == g
    Z[ix, ] <- abs(sweep(R[ix, , drop = FALSE], 2, col_med(R[ix, , drop = FALSE])))
  }
  Z
}

# Brown-Forsythe with two groups is a two-sample t on the spread variable
bf_p <- function(Z) {
  a <- Z[is_t, , drop = FALSE]; b <- Z[!is_t, , drop = FALSE]
  na <- nrow(a); nb <- nrow(b)
  ma <- colMeans(a); mb <- colMeans(b)
  va <- colSums(sweep(a, 2, ma)^2) / (na - 1)
  vb <- colSums(sweep(b, 2, mb)^2) / (nb - 1)
  sp <- ((na - 1) * va + (nb - 1) * vb) / (na + nb - 2)
  t  <- (ma - mb) / sqrt(sp * (1/na + 1/nb))
  2 * pt(abs(t), df = na + nb - 2, lower.tail = FALSE)
}

disp_p <- function(Z) {
  f  <- lm.fit(Xdis, Z)
  s2 <- colSums(f$residuals^2) / (n - ncol(Xdis))
  XtXi <- chol2inv(chol(crossprod(Xdis)))
  t  <- f$coefficients[jdis, ] / sqrt(s2 * XtXi[jdis, jdis])
  2 * pt(abs(t), df = n - ncol(Xdis), lower.tail = FALSE)
}

simulate <- function(sd_vec, label) {
  # sd_vec: per-donor residual SD. Y is pure noise with that SD structure,
  # so there is no mean signal and no gene is "real".
  Y <- matrix(rnorm(n * NGENE), n, NGENE) * sd_vec
  R <- lm.fit(Xloc, Y)$residuals      # same location model as the manuscript
  Z <- spread(R)
  p_bf <- bf_p(Z); p_dm <- disp_p(Z)
  c(BF = mean(p_bf < ALPHA, na.rm = TRUE),
    DISP = mean(p_dm < ALPHA, na.rm = TRUE))
}

cat("design: n =", n, " T2D =", sum(is_t), " ND =", sum(!is_t), "\n")
cat("correlation of beta score with group:",
    round(cor(zb, as.numeric(is_t)), 3), "\n")
cat("correlation of imm  score with group:",
    round(cor(zi, as.numeric(is_t)), 3), "\n\n")

## Scenario A -- flat
A <- simulate(rep(1, n), "A")

## Scenario B -- dispersion depends on composition, NOT on group.
## sd multiplies by exp(gamma * z_beta); gamma = 0.35 gives roughly the
## spread heterogeneity seen in the real residuals. Two strengths shown.
B1 <- simulate(exp(-0.35 * zb), "B1")
B2 <- simulate(exp(-0.50 * zb + 0.20 * zi), "B2")

## Scenario C -- a genuine group effect of the size claimed (SD ratio 2)
C <- simulate(ifelse(is_t, 2, 1), "C")

tab <- rbind("A  null, no composition dependence"      = A,
             "B1 null, dispersion ~ beta (moderate)"   = B1,
             "B2 null, dispersion ~ beta + imm"        = B2,
             "C  TRUE group effect, SD ratio 2"        = C)
colnames(tab) <- c("Brown-Forsythe", "dispersion model")

cat("===== rejection rate at alpha = 0.05 =====\n")
print(round(tab, 4))

cat("\n===== HOW TO READ THIS =====\n")
cat("  Rows A, B1, B2 contain NO group effect. Any rate above ~0.05 there\n")
cat("  is a false positive rate, not power.\n\n")
cat("  If Brown-Forsythe is inflated in B1/B2 and the dispersion model is\n")
cat("  near 0.05, then composition alone can manufacture the BF result and\n")
cat("  the published UNC13D P = 0.028 is not trustworthy. The dispersion\n")
cat("  model becomes the test to report.\n\n")
cat("  If both hold ~0.05 in B1/B2, the two tests answer different questions\n")
cat("  rather than one being wrong, and both belong in the paper with the\n")
cat("  distinction stated.\n\n")
cat("  If the dispersion model is BELOW 0.05 in A and has low power in C,\n")
cat("  it is conservative and its null result for UNC13D is weaker evidence\n")
cat("  than it looks. Row C is what tells you that.\n")

saveRDS(tab, "C_calibration.rds")
write.csv(tab, "C_calibration.csv")
cat("\nWrote C_calibration.csv\n")
