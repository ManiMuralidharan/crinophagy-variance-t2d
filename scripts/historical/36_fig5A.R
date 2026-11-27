# =====================================================================
# 36. Rebuild Figure 5A (cross-platform replication)
#
# ONLY CHANGE: the discovery axis label, n = 123 -> n = 103.
# Per-protein analyses used the 103 donors in the two-group comparison
# (87 non-diabetic, 16 T2D), per Methods 4.2 and the S6 legend.
#
# The replication estimates are recomputed from session_full.RData.
# The discovery estimates CANNOT be recomputed here, because that
# session holds the replication cohort only. Fill them in below from
# Table S1 or from the earlier discovery output. Do not guess them.
#
# RUN:  source("36_fig5A.R")
# =====================================================================

library(ggplot2); library(dplyr); library(car)

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))
if (!exists("mat")) load("session_full.RData")
stopifnot(exists("mat"), exists("pd2"))
set.seed(1)


# ---- 1. DISCOVERY VALUES: fill these in -----------------------------
# From Table S1 / the discovery per-protein output. The point estimates
# below are the published ones. REPLACE the CI values with the real
# bootstrap limits before using the figure.

discovery <- data.frame(
  gene  = c("UNC13D", "RAB3A"),
  ratio = c(2.065, 1.725),
  ci_lo = c(NA, NA),        # <- FILL IN
  ci_hi = c(NA, NA),        # <- FILL IN
  cohort = "Discovery\n(proteomics, n = 103)"
)

if (any(is.na(discovery$ci_lo))) {
  cat("\n*** Discovery confidence intervals are not filled in.\n")
  cat("*** Find them in Table S1 (columns ci_lo / ci_hi for UNC13D\n")
  cat("*** and RAB3A in the discovery cohort) and edit section 1.\n")
  cat("*** Stopping so the figure is not drawn with missing bars.\n\n")
  stop("discovery CIs required")
}


# ---- 2. REPLICATION VALUES: recomputed ------------------------------

TARGETS <- c(UNC13D = "ENSG00000092929", RAB3A = "ENSG00000105649")
stopifnot(all(TARGETS %in% rownames(mat)))

one <- function(ens, nboot = 5000) {
  d <- pd2
  d$y <- mat[ens, rownames(d)]
  d$r <- residuals(lm(y ~ batch + imm + beta, data = d))
  a <- d$r[d$group == "T2D"]; b <- d$r[d$group == "ND"]
  bt <- replicate(nboot,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))
  c(ratio = sd(a) / sd(b),
    ci_lo = unname(quantile(bt, .025)),
    ci_hi = unname(quantile(bt, .975)))
}

rep_vals <- t(sapply(TARGETS, one))
replication <- data.frame(
  gene   = names(TARGETS),
  ratio  = rep_vals[, "ratio"],
  ci_lo  = rep_vals[, "ci_lo"],
  ci_hi  = rep_vals[, "ci_hi"],
  cohort = sprintf("Replication\n(RNA-seq, n = %d)", nrow(pd2)),
  row.names = NULL
)

cat("\n===== RECOMPUTED REPLICATION VALUES =====\n")
print(replication, digits = 4, row.names = FALSE)
cat("\nPublished values are UNC13D 2.001 and RAB3A 1.101.\n")
cat("If these differ, STOP and reconcile before redrawing.\n\n")


# ---- 3. Draw --------------------------------------------------------

D <- bind_rows(discovery, replication)
D$cohort <- factor(D$cohort, levels = c(unique(discovery$cohort),
                                        unique(replication$cohort)))
D$gene <- factor(D$gene, levels = c("RAB3A", "UNC13D"))

p <- ggplot(D, aes(cohort, ratio, colour = gene, group = gene)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_line(linewidth = .6, position = position_dodge(width = .12)) +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi), size = .5,
                  linewidth = .7, position = position_dodge(width = .12)) +
  scale_colour_manual(values = c(UNC13D = "#1b5e20", RAB3A = "#c62828"),
                      name = NULL) +
  scale_y_continuous(limits = c(0.9, 3.0)) +
  labs(x = NULL, y = "SD ratio (T2D / non-diabetic), 95% CI",
       title = "A   Cross-platform replication") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(.88, .88),
        plot.title = element_text(size = 10, face = "bold"))

ggsave("figures/Fig5A_corrected.pdf", p, width = 3.4, height = 3.4)
cat("Wrote figures/Fig5A_corrected.pdf\n")
cat("Check it against the published panel before replacing it.\n")
