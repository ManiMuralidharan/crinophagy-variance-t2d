# =====================================================================
# 41_between_set.R   C-1 and C-2
#
# Answers, with one procedure:
#   R1-M2  "the random-set test does not test disease association"
#   R2-M8  "a difference of P values is not a test of a difference"
#
# Design: permute glycemic labels WITHIN BATCH and recompute every SD
# ratio. Because no covariate model contains `group`, the residuals are
# unchanged by the permutation, so this costs one matrix operation per
# replicate. The null preserves gene-gene correlation, expression level,
# detectability and baseline variability exactly, because it is the same
# genes every time. That is what the size-matched gene-resampling null
# could not do.
#
# RUN:  source("40_common.R"); source("41_between_set.R")
# Time: a few minutes at NPERM = 1e5.
# =====================================================================

NPERM <- 1e5          # raise to 1e6 for the final numbers
set.seed(1)

stopifnot(exists("RES2"), exists("permute_labels"))

crino_i <- match(crino, colnames(RES2))
macro_i <- match(macro, colnames(RES2))
stopifnot(!anyNA(crino_i), !anyNA(macro_i))
cat(sprintf("Crinophagy genes: %d   Comparator genes: %d\n",
            length(crino_i), length(macro_i)))

run_model <- function(RES, label) {

  r_obs  <- sd_ratios(RES, OBS_LAB)
  m_cri  <- median(r_obs[crino_i])
  m_mac  <- median(r_obs[macro_i])
  d_obs  <- m_cri - m_mac

  null_cri <- numeric(NPERM)
  null_mac <- numeric(NPERM)
  null_dif <- numeric(NPERM)

  pb <- txtProgressBar(min = 0, max = NPERM, style = 3)
  for (i in seq_len(NPERM)) {
    r <- sd_ratios(RES, permute_labels())
    null_cri[i] <- median(r[crino_i])
    null_mac[i] <- median(r[macro_i])
    null_dif[i] <- null_cri[i] - null_mac[i]
    if (i %% 2000 == 0) setTxtProgressBar(pb, i)
  }
  close(pb)

  # one-sided for the set tests (the hypothesis is variance EXPANSION),
  # two-sided for the difference (we did not prespecify its direction)
  p_cri <- (sum(null_cri >= m_cri) + 1) / (NPERM + 1)
  p_mac <- (sum(null_mac >= m_mac) + 1) / (NPERM + 1)
  p_dif <- (sum(abs(null_dif) >= abs(d_obs)) + 1) / (NPERM + 1)

  cat(sprintf("\n===== %s =====\n", label))
  cat(sprintf("  crinophagy  median SD ratio %.4f   null median %.4f   label-permutation P %.5f\n",
              m_cri, median(null_cri), p_cri))
  cat(sprintf("  comparator  median SD ratio %.4f   null median %.4f   label-permutation P %.5f\n",
              m_mac, median(null_mac), p_mac))
  cat(sprintf("  DIFFERENCE  %.4f  (95%% null range %.4f to %.4f)   two-sided P %.5f\n",
              d_obs, quantile(null_dif, .025), quantile(null_dif, .975), p_dif))
  cat("\n  Note: these P values test DISEASE ASSOCIATION (labels permuted),\n")
  cat("  which is a different question from the published gene-resampling\n")
  cat("  P values (0.0083 / 0.036), which test extremeness among random sets\n")
  cat("  of the same size. Report both and say which is which.\n")

  list(label = label, m_cri = m_cri, m_mac = m_mac, d_obs = d_obs,
       p_cri = p_cri, p_mac = p_mac, p_dif = p_dif,
       null_cri = null_cri, null_mac = null_mac, null_dif = null_dif)
}

r_m1 <- run_model(RES1, "Model 1, batch only")
r_m2 <- run_model(RES2, "Model 2, fully adjusted (primary)")

save(r_m1, r_m2, NPERM, file = "C1_C2_between_set.RData")

out <- data.frame(
  model          = c(r_m1$label, r_m2$label),
  median_crino   = c(r_m1$m_cri,  r_m2$m_cri),
  median_comp    = c(r_m1$m_mac,  r_m2$m_mac),
  difference     = c(r_m1$d_obs,  r_m2$d_obs),
  P_crino_label  = c(r_m1$p_cri,  r_m2$p_cri),
  P_comp_label   = c(r_m1$p_mac,  r_m2$p_mac),
  P_difference   = c(r_m1$p_dif,  r_m2$p_dif)
)
write.csv(out, "C1_C2_between_set.csv", row.names = FALSE)
cat("\nWrote C1_C2_between_set.csv and .RData\n")

cat("\nINTERPRET HONESTLY:\n")
cat("  If P_difference is NOT significant, you cannot claim the crinophagy set\n")
cat("  differs from the comparator set. You can still claim the crinophagy set\n")
cat("  is extreme genome-wide and that the comparator is not, which is a weaker\n")
cat("  but defensible statement. Write whichever one the number supports.\n")
