# =====================================================================
# 57_pathway_ci.R      RESOLVES R2-m4
#
# "Explain why the discovery pathway SD-ratio confidence interval
#  includes 1 (0.97-2.21) despite the Brown-Forsythe P value of 0.018."
#
# This is about the DISCOVERY ssGSEA crinophagy score, not about any
# gene. The two numbers are not in conflict; they estimate different
# things. This script reproduces both from your own data and quantifies
# why they differ, so the answer in the response letter is demonstrated
# rather than asserted.
#
# RUN:  source("57_pathway_ci.R")     (does not need 40_common.R)
# =====================================================================

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))
set.seed(1)

PUB <- list(ratio = 1.61, ci = c(0.97, 2.21), bf_p = 0.018)
NBOOT <- 20000

ss  <- read.csv("ssgsea_scores_with_metadata.csv", stringsAsFactors = FALSE)
don <- read.csv("donor.csv", stringsAsFactors = FALSE)
cat(sprintf("ssGSEA scores: %d donors\n", nrow(ss)))

# islet recovery is not in the ssGSEA export; bring it from donor.csv
ss$recovery <- don$percentieqrecoverypostculture[match(ss$record_id, don$record_id)]
ss$dx_comp  <- don$diagnosis_computed[match(ss$record_id, don$record_id)]
cat("diagnosis_computed in the score file:\n"); print(table(ss$dx_comp, useNA = "ifany"))

ss$grp <- ifelse(ss$dx_comp == "Type2", "T2D",
          ifelse(ss$dx_comp == "Pre.T2D", "preT2D", "ND"))

COV <- c("age", "sex", "bmi", "culture_time", "recovery")
ok  <- complete.cases(ss[, c(COV, "crinophagy")]) & ss$grp %in% c("ND", "T2D")
d   <- ss[ok, ]
cat(sprintf("\ncovariate-complete two-group donors: %d\n", nrow(d)))
print(table(d$grp))
cat("published: 87 non-diabetic, 16 T2D\n")

d$sex <- factor(d$sex)
fit <- lm(crinophagy ~ age + sex + bmi + culture_time + recovery, data = d)
d$r <- residuals(fit)

a <- d$r[d$grp == "T2D"]; b <- d$r[d$grp == "ND"]
ratio <- sd(a) / sd(b)

bt <- replicate(NBOOT, sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))
ci <- quantile(bt, c(.025, .975))

z  <- abs(d$r - ave(d$r, d$grp, FUN = median))
bf <- t.test(z[d$grp == "T2D"], z[d$grp == "ND"], var.equal = TRUE)

cat("\n===== SELF-CHECK against the published pathway numbers =====\n")
cat(sprintf("  SD ratio      reconstructed %.3f    published %.2f\n", ratio, PUB$ratio))
cat(sprintf("  95%% CI        reconstructed [%.2f, %.2f]  published [%.2f, %.2f]\n",
            ci[1], ci[2], PUB$ci[1], PUB$ci[2]))
cat(sprintf("  Brown-Forsythe reconstructed P = %.4f   published %.3f\n",
            bf$p.value, PUB$bf_p))
okchk <- abs(ratio - PUB$ratio) < 0.08 && abs(bf$p.value - PUB$bf_p) < 0.01
cat(ifelse(okchk, "  Self-check passed.\n",
           "  *** MISMATCH. Check the covariate set and the score column before quoting anything.\n"))

cat("\n===== WHY THE TWO DISAGREE =====\n")
cat(sprintf("  1. Different estimands.\n"))
cat(sprintf("     Interval is for a RATIO of SDs            : %.3f  [%.2f, %.2f]\n",
            ratio, ci[1], ci[2]))
cat(sprintf("     Brown-Forsythe tests a DIFFERENCE in mean |deviation|:\n"))
cat(sprintf("       T2D %.4f vs ND %.4f, difference CI [%.4f, %.4f], P = %.4f\n",
            mean(z[d$grp == "T2D"]), mean(z[d$grp == "ND"]),
            -bf$conf.int[2], -bf$conf.int[1], bf$p.value))
cat(sprintf("\n  2. The ratio's sampling distribution is right-skewed at n = %d vs %d.\n",
            length(a), length(b)))
cat(sprintf("     bootstrap skewness %.2f; median %.3f against point estimate %.3f\n",
            mean((bt - mean(bt))^3) / sd(bt)^3, median(bt), ratio))
cat(sprintf("     percentile interval is wide on the low side by construction.\n"))
cat(sprintf("\n  3. The bootstrap resamples donors; Brown-Forsythe conditions on them.\n"))
cat(sprintf("     With only %d T2D donors the bootstrap carries resampling variability\n", length(a)))
cat(sprintf("     that the fixed-sample test does not.\n"))

# how much of the interval width is driven by the small T2D group
bt_nd <- replicate(5000, sd(a) / sd(sample(b, length(b), TRUE)))
bt_t2 <- replicate(5000, sd(sample(a, length(a), TRUE)) / sd(b))
cat(sprintf("\n  width with only ND resampled  : %.2f\n", diff(quantile(bt_nd, c(.025, .975)))))
cat(sprintf("  width with only T2D resampled : %.2f\n", diff(quantile(bt_t2, c(.025, .975)))))
cat(sprintf("  width with both               : %.2f\n", diff(ci)))
cat("  The larger of the first two is where the imprecision comes from.\n")

## leave-one-out, to support the Figure S1 statement
loo <- sapply(which(d$grp == "T2D"), function(i) {
  dd <- d[-i, ]
  dd$r <- residuals(lm(crinophagy ~ age + sex + bmi + culture_time + recovery, data = dd))
  zz <- abs(dd$r - ave(dd$r, dd$grp, FUN = median))
  t.test(zz[dd$grp == "T2D"], zz[dd$grp == "ND"], var.equal = TRUE)$p.value
})
names(loo) <- d$record_id[d$grp == "T2D"]
cat("\n===== LEAVE-ONE-OUT over the 16 T2D donors (supports Figure S1) =====\n")
cat(sprintf("  full-data P = %.4f; range after dropping one donor: %.4f to %.4f\n",
            bf$p.value, min(loo), max(loo)))
cat(sprintf("  donors whose removal pushes P above 0.05: %d of %d\n",
            sum(loo > 0.05), length(loo)))
print(round(sort(loo, decreasing = TRUE)[1:5], 4))

out <- data.frame(statistic = c("SD ratio", "CI low", "CI high", "Brown-Forsythe P",
                                "max leave-one-out P", "n above 0.05"),
                  value = c(ratio, ci[1], ci[2], bf$p.value, max(loo), sum(loo > 0.05)))
write.csv(out, "C_pathway_ci_vs_bf.csv", row.names = FALSE)
cat("\nWrote C_pathway_ci_vs_bf.csv\n")

cat("\n===== DRAFT SENTENCE FOR THE RESPONSE =====\n")
cat("  The interval and the test target different quantities. The bootstrap\n")
cat("  interval is for the ratio of standard deviations, whose sampling\n")
cat("  distribution is right-skewed at these group sizes, whereas the\n")
cat("  Brown-Forsythe statistic tests a difference in mean absolute deviation\n")
cat("  from the group median and conditions on the observed donors. We now\n")
cat("  report both with their estimands stated and do not present either as\n")
cat("  confirming the other.\n")
