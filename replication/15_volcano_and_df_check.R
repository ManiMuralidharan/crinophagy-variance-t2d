# =====================================================================
# Volcano plot + degrees-of-freedom sensitivity check
#
# TWO SEPARATE THINGS:
#
# (A) DF CORRECTION - this is a real methodological issue raised in
#     review and must be settled before the manuscript is final.
#     Levene's test and the SD ratio use n-1 in the variance
#     denominator, but Model 2 estimates p parameters (intercept,
#     batch levels, immune score, beta-cell score). Residual variance
#     is therefore downward-biased, and because n differs between
#     groups (18 ND vs 36 T2D) the bias is unequal. This recomputes
#     every ratio with an n-p denominator and reports whether the
#     ordering or the conclusions change.
#
#     Expectation: a small uniform shift, no change in ranking. If the
#     ranking DOES change, the manuscript numbers must be revised to
#     the corrected version, not the other way round.
#
# (B) VOLCANO - supplementary figure, modelled on panel D of the
#     uploaded reference (focused pathway plot, colour-coded by
#     functional subgroup) rather than panel B (genome-wide scatter).
#     Note that with 46 genes and few significant results this plot is
#     sparse; Fig 4 (forest) carries the same information more
#     honestly. Supplementary, not main.
#
# Requires: mat, pd2, bg2 in session; canonical_full.csv on disk.
# =====================================================================

library(ggplot2)
library(dplyr)
library(car)
library(ggrepel)

setwd("/Users/mmani3/Downloads/Human Islets")
REP  <- "replication_GSE164416"
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

full <- read.csv(file.path(REP, "canonical_full.csv"),
                 stringsAsFactors = FALSE)
m2 <- full %>% dplyr::filter(model == "M2_full")


# =====================================================================
# (A) DEGREES-OF-FREEDOM SENSITIVITY
# =====================================================================

cat("===== DF CORRECTION SENSITIVITY =====\n\n")

# Model 2 design matrix: how many parameters are actually estimated?
X <- model.matrix(~ batch + imm + beta, data = pd2)
p_est <- qr(X)$rank
n_tot <- nrow(pd2)
n_t2d <- sum(pd2$group == "T2D")
n_nd  <- sum(pd2$group == "ND")

cat(sprintf("Donors: %d (ND %d, T2D %d)\n", n_tot, n_nd, n_t2d))
cat(sprintf("Parameters estimated in Model 2: %d\n", p_est))
cat(sprintf("Residual df: %d\n\n", n_tot - p_est))

# Allocate the df loss to each group in proportion to its size. This is
# an approximation: residual df are not cleanly separable by group, but
# it is the conservative reading of the reviewer's point and shows
# whether the concern is material.
df_loss_t2d <- (p_est - 1) * n_t2d / n_tot
df_loss_nd  <- (p_est - 1) * n_nd  / n_tot

cat(sprintf("Approx. df allocated: ND %.2f, T2D %.2f\n\n",
            df_loss_nd, df_loss_t2d))

df_check <- function(ens, label) {
  d <- pd2
  d$y <- mat[ens, rownames(d)]
  d$r <- residuals(lm(y ~ batch + imm + beta, data = d))

  a <- d$r[d$group == "T2D"]
  b <- d$r[d$group == "ND"]

  # standard: n-1 denominator
  v_a1 <- sum(a^2) / (length(a) - 1)
  v_b1 <- sum(b^2) / (length(b) - 1)

  # df-corrected: subtract the group's share of estimated parameters
  v_a2 <- sum(a^2) / (length(a) - 1 - df_loss_t2d)
  v_b2 <- sum(b^2) / (length(b) - 1 - df_loss_nd)

  data.frame(gene = label,
             ratio_n1 = sqrt(v_a1 / v_b1),
             ratio_np = sqrt(v_a2 / v_b2),
             row.names = NULL)
}

key <- m2 %>% dplyr::select(gene, ensembl) %>% dplyr::filter(!is.na(ensembl))
dfres <- do.call(rbind, lapply(seq_len(nrow(key)), function(i)
  df_check(key$ensembl[i], key$gene[i])))

dfres <- dfres %>%
  mutate(diff = ratio_np - ratio_n1,
         pct_change = 100 * diff / ratio_n1) %>%
  arrange(desc(ratio_np))

cat("Effect of the correction on every panel gene:\n")
print(dfres, digits = 4, row.names = FALSE)

cat(sprintf("\nMean absolute change: %.4f (%.2f%%)\n",
            mean(abs(dfres$diff)), mean(abs(dfres$pct_change))))

# Does the ranking change?
r1 <- rank(-dfres$ratio_n1); r2 <- rank(-dfres$ratio_np)
cat(sprintf("Spearman correlation of rankings: %.6f\n",
            cor(r1, r2, method = "spearman")))
cat(sprintf("Genes changing rank: %d of %d\n",
            sum(r1 != r2), nrow(dfres)))

top <- dfres %>% dplyr::filter(gene %in% c("UNC13D","STX7","CTSD","GABARAP"))
cat("\nKey genes:\n"); print(top, digits = 4, row.names = FALSE)

cat("\nINTERPRETATION\n")
if (cor(r1, r2, method = "spearman") > 0.999 &&
    mean(abs(dfres$pct_change)) < 2) {
  cat("  The correction applies an essentially uniform scaling. Ranking\n")
  cat("  and percentiles are unaffected, because the ratio is a quotient\n")
  cat("  of two variances estimated from the same model. Report the\n")
  cat("  uncorrected values and state that df correction was checked\n")
  cat("  and changed nothing material.\n")
} else {
  cat("  The correction is NOT negligible. Revise the manuscript to the\n")
  cat("  df-corrected values throughout, and report both.\n")
}

write.csv(dfres, file.path(REP, "df_sensitivity.csv"), row.names = FALSE)


# =====================================================================
# (B) VOLCANO PLOT (supplementary)
# =====================================================================

# Functional subgroups, for colour only. These are descriptive
# annotations of the frozen sets, NOT a tested partition: the
# fusion/cargo distinction was suggested by the results and is not
# part of the pre-registered design. Colour communicates biology to
# the reader; no statistic is computed on these groups.

grp <- c(
  UNC13D="Fusion", STX7="Fusion", STX8="Fusion", STX17="Fusion",
  VAMP7="Fusion", VAMP8="Fusion", SNAP23="Fusion", VPS33A="Fusion",
  VPS16="Fusion", VPS11="Fusion",
  CTSB="Lysosomal", CTSD="Lysosomal", CTSL="Lysosomal",
  LAMP1="Lysosomal", LAMP2="Lysosomal", ATP6V0D1="Lysosomal",
  ATP6V1A="Lysosomal", TFEB="Lysosomal", PIKFYVE="Lysosomal",
  INS="Cargo", CHGA="Cargo", CHGB="Cargo", PCSK1="Cargo",
  PCSK2="Cargo", CPE="Cargo",
  RAB3A="Granule", RAB3GAP1="Granule", RAB27A="Granule",
  CADPS="Granule", SYTL4="Granule", PCLO="Granule"
)

v <- m2 %>%
  mutate(
    subgroup = ifelse(set == "comparator", "Macroautophagy",
               ifelse(set == "post-hoc", "Post-hoc comparator",
               ifelse(gene %in% names(grp), grp[gene], "Other"))),
    log2ratio = log2(sd_ratio),
    negl10p = -log10(p),
    lab = gene %in% c("UNC13D","STX7","CTSD","GABARAP","RAB3A",
                      "INS","STX1A","RAB27A","MAP1LC3B","VPS16")
  )

thr99 <- log2(quantile(bg2, .99))

pal <- c(Fusion = "#1b5e20", Lysosomal = "#2e7d32", Cargo = "#ef6c00",
         Granule = "#8d6e63", Macroautophagy = "#90a4ae",
         `Post-hoc comparator` = "#c62828", Other = "grey60")

vol <- ggplot(v, aes(log2ratio, negl10p)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf,
           ymax = -log10(0.05), fill = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = thr99, linetype = "dotted", colour = "grey20") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted",
             colour = "grey20") +
  geom_point(aes(fill = subgroup), shape = 21, size = 2.6,
             colour = "grey20", stroke = .3) +
  geom_text_repel(data = dplyr::filter(v, lab), aes(label = gene),
                  size = 2.7, min.segment.length = 0,
                  box.padding = .45, seed = 1, max.overlaps = 20) +
  scale_fill_manual(values = pal, name = NULL) +
  annotate("text", x = -0.35, y = 0.06, label = "lower variance in T2D",
           size = 2.6, colour = "grey45", hjust = 0) +
  annotate("text", x = 0.35, y = 0.06, label = "higher variance in T2D",
           size = 2.6, colour = "grey45", hjust = 0) +
  labs(x = expression(paste("Effect size (", log[2],
                            " variance ratio, T2D / non-diabetic)")),
       y = expression(paste(-log[10], "(Levene P)")),
       title = "Variance heterogeneity across the pre-registered panel",
       subtitle = "Fully adjusted model. Dotted vertical line: 99th percentile of the genome-wide background",
       caption = paste("Colours are descriptive functional annotations of the frozen gene sets, not a tested",
                       "\npartition. P values are uncorrected; Benjamini-Hochberg correction was applied within",
                       "\nthe comparator panel only (see Methods). Grey band: P > 0.05.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        plot.caption = element_text(size = 7, colour = "grey40", hjust = 0),
        legend.position = "right")

ggsave(file.path(FIGS, "FigS_volcano_variance.pdf"), vol,
       width = 7.5, height = 5.5)

cat("\nWrote FigS_volcano_variance.pdf and df_sensitivity.csv\n")
cat("\nNote: this volcano is sparse by nature - few genes reach P<0.05\n")
cat("at n = 18 vs 36. Fig 4 (forest with CIs) conveys the same result\n")
cat("more honestly. Keep the volcano supplementary.\n")
