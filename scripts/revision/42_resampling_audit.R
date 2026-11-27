# =====================================================================
# 42_resampling_audit.R   C-7 and C-6
#
# C-7  Does the bootstrap refit the nuisance model inside each replicate?
#      If it does not, the intervals are too narrow and Referee 2 is
#      entitled to say so. This script quantifies the difference rather
#      than asserting there is none.
#
# C-6  Coverage of the bootstrap interval at the observed group sizes.
#
# RUN:  source("40_common.R"); source("42_resampling_audit.R")
# =====================================================================

NBOOT <- 5000
NSIM  <- 2000         # coverage replicates; raise if time allows
set.seed(1)

TARGETS <- c(UNC13D = "ENSG00000092929",
             STX7   = "ENSG00000079950",
             CTSD   = "ENSG00001xxxxx")        # <-- FILL IN CTSD's Ensembl ID
TARGETS <- TARGETS[TARGETS %in% rownames(mat)]
if (!length(TARGETS)) stop("No target gene IDs matched rownames(mat). Fill TARGETS in.")
cat("Targets resolved:", paste(names(TARGETS), collapse = ", "), "\n")


## ---- C-7: bootstrap WITHOUT vs WITH nuisance refit ------------------

boot_norefit <- function(ens, nboot = NBOOT) {
  # residuals fitted once on the full data, then resampled (the common,
  # and optimistic, implementation)
  r <- RES2[, ens]
  a <- r[OBS_LAB == "T2D"]; b <- r[OBS_LAB == "ND"]
  replicate(nboot,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))
}

boot_refit <- function(ens, nboot = NBOOT) {
  # donors resampled FIRST, nuisance model refitted inside each replicate
  y   <- mat[ens, SAMPLES]
  idx <- seq_len(nrow(pd2))
  out <- numeric(nboot)
  for (i in seq_len(nboot)) {
    s <- sample(idx, length(idx), TRUE)
    d <- pd2[s, , drop = FALSE]; d$y <- y[s]
    # a resample can lose a batch level or a group; skip degenerate draws
    if (length(unique(d$group)) < 2 || any(table(d$group) < 3)) { out[i] <- NA; next }
    d$batch <- droplevels(factor(d$batch))
    fit <- try(lm(y ~ batch + imm + beta, data = d), silent = TRUE)
    if (inherits(fit, "try-error")) { out[i] <- NA; next }
    r  <- residuals(fit)
    out[i] <- sd(r[d$group == "T2D"]) / sd(r[d$group == "ND"])
  }
  out
}

cat("\n===== C-7  effect of refitting the nuisance model =====\n")
c7 <- do.call(rbind, lapply(names(TARGETS), function(g) {
  e  <- TARGETS[[g]]
  bn <- boot_norefit(e); br <- boot_refit(e)
  qn <- quantile(bn, c(.025, .975), na.rm = TRUE)
  qr <- quantile(br, c(.025, .975), na.rm = TRUE)
  data.frame(gene = g,
             point      = sd_ratios(RES2, OBS_LAB)[e],
             lo_norefit = qn[1], hi_norefit = qn[2], width_norefit = diff(qn),
             lo_refit   = qr[1], hi_refit   = qr[2], width_refit   = diff(qr),
             pct_wider  = 100 * (diff(qr) / diff(qn) - 1),
             failed_draws = sum(is.na(br)))
}))
print(c7, row.names = FALSE, digits = 4)
cat("\n  If pct_wider is small (say under 10%), say in the Methods that the\n")
cat("  published intervals do not refit and that refitting changes them\n")
cat("  negligibly, and give this table. If it is large, use the refitted\n")
cat("  intervals in the manuscript. Do not report only the narrower one.\n")


## ---- C-6: coverage at the observed design ---------------------------
# Simulate from the fitted residual distribution with a KNOWN SD ratio
# and ask how often the nominal 95% interval contains it.

cat("\n===== C-6  bootstrap interval coverage at the observed group sizes =====\n")
n_t2d <- sum(OBS_LAB == "T2D"); n_nd <- sum(OBS_LAB == "ND")
cat(sprintf("  group sizes: %d T2D, %d non-diabetic\n", n_t2d, n_nd))

coverage <- function(true_ratio, nsim = NSIM, nboot = 1000) {
  hit <- logical(nsim)
  for (i in seq_len(nsim)) {
    a <- rnorm(n_t2d, 0, true_ratio); b <- rnorm(n_nd, 0, 1)
    bt <- replicate(nboot,
      sd(sample(a, n_t2d, TRUE)) / sd(sample(b, n_nd, TRUE)))
    q <- quantile(bt, c(.025, .975))
    hit[i] <- q[1] <= true_ratio && true_ratio <= q[2]
  }
  mean(hit)
}

c6 <- data.frame(true_sd_ratio = c(1.0, 1.5, 2.0))
c6$coverage <- sapply(c6$true_sd_ratio, coverage)
print(c6, row.names = FALSE, digits = 3)
cat("\n  Nominal is 0.95. Percentile bootstrap at these group sizes is\n")
cat("  usually a little below nominal. Report the number you get, whatever\n")
cat("  it is; an honest coverage figure is a stronger answer to R2-M2 than\n")
cat("  an assertion that the method is adequate.\n")

save(c7, c6, file = "C6_C7_resampling_audit.RData")
write.csv(c7, "C7_refit_comparison.csv", row.names = FALSE)
write.csv(c6, "C6_coverage.csv", row.names = FALSE)
cat("\nWrote C6_coverage.csv, C7_refit_comparison.csv, C6_C7_resampling_audit.RData\n")
