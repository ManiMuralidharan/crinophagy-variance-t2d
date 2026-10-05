# =====================================================================
# PUBLICATION FIGURES - five PDFs from frozen outputs only
#
# This script performs NO new statistics. Every value plotted is read
# from a frozen file or from bg1/bg2 computed in script 13. If a number
# here disagrees with the manuscript, the manuscript is wrong, not this.
#
#   Fig 1  Genome-wide variance distribution (the empirical null)
#   Fig 2  Model 1 -> Model 2 dumbbell (composition adjustment)
#   Fig 3  Mean vs variance quadrant map (cross-cohort; see caption)
#   Fig 4  Forest plot, all 46 genes, fully adjusted
#   Fig 5  Cross-platform replication: UNC13D vs RAB3A
#
# DELIBERATELY NOT BUILT: gene-gene correlation matrices and ROC
# classifier panels. Correlation estimates at n = 18 vs 36 are far too
# unstable to support a "network rewiring" claim, and a classifier
# panel invites a diagnostic-utility claim this sample size cannot
# carry. Both also sit outside the frozen analysis.
#
# Requires: session_gse164416.RData (mat, pd2, ratio_all) and, if
# regenerating the two backgrounds, script 13 sections 2-3.
# =====================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

setwd("/Users/mmani3/Downloads/Human Islets")
REP  <- "replication_GSE164416"
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

# consistent look across all five
theme_pub <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey30"),
        plot.caption = element_text(size = 7, colour = "grey40", hjust = 0))

COL <- c(crinophagy = "#1b5e20", comparator = "#90a4ae",
         `post-hoc` = "#ef6c00")

cmp <- read.csv(file.path(REP, "canonical_both_models.csv"),
                stringsAsFactors = FALSE)
full <- read.csv(file.path(REP, "canonical_full.csv"),
                 stringsAsFactors = FALSE)
m2 <- full %>% filter(model == "M2_full")

if (!exists("bg2")) {
  stop("bg1/bg2 not in session. Re-run script 13 sections 2-3 first, ",
       "or load a saved copy.")
}


# ---------------------------------------------------------------------
# FIG 1 - genome-wide variance distribution
# Answers the first question any reviewer asks: isn't everything more
# variable in T2D? The distribution centres near 1.0; the marked genes
# are in the tail.
# ---------------------------------------------------------------------

d1 <- data.frame(ratio = as.numeric(bg2))
mark <- m2 %>%
  filter(gene %in% c("UNC13D","STX7","CTSD","GABARAP","RAB3A","INS")) %>%
  select(gene, sd_ratio, set)

p99 <- quantile(bg2, .99)

f1 <- ggplot(d1, aes(ratio)) +
  geom_histogram(bins = 120, fill = "grey85", colour = NA) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_vline(xintercept = p99, linetype = "dotted", colour = "grey20") +
  geom_vline(data = mark, aes(xintercept = sd_ratio, colour = set),
             linewidth = .5) +
  geom_text(data = mark,
            aes(x = sd_ratio, y = Inf, label = gene, colour = set),
            angle = 90, hjust = 1.1, vjust = -0.4, size = 2.6) +
  scale_colour_manual(values = COL, guide = "none") +
  coord_cartesian(xlim = c(0.3, 3)) +
  labs(x = "SD ratio (T2D / non-diabetic)", y = "Genes",
       title = "Genome-wide distribution of inter-donor variance ratios",
       subtitle = sprintf("%s expressed genes, fully adjusted model. Dotted line: 99th percentile (%.2f)",
                          format(length(bg2), big.mark = ","), p99),
       caption = "Adjusted for library preparation batch, immune infiltration score and beta-cell content score.") +
  theme_pub

ggsave(file.path(FIGS, "Fig1_genomewide_distribution.pdf"), f1,
       width = 6.5, height = 4)


# ---------------------------------------------------------------------
# FIG 2 - composition adjustment
# The argument: adjustment removes composition-driven variance where it
# exists (GABARAP, STX1A, MAP1LC3A, INS attenuate) and does not remove
# the crinophagy signals (UNC13D, STX7, CTSD hold or rise).
# ---------------------------------------------------------------------

sel <- c("UNC13D","STX7","CTSD","CTSB","GABARAP","GABARAPL2",
         "STX1A","MAP1LC3A","INS","RAB3A","SNAP23","CPE",
         "ATP6V1A","CALCOCO2","ATG3","RAB27A")

d2 <- cmp %>%
  filter(gene %in% sel) %>%
  mutate(gene = reorder(gene, r2),
         direction = ifelse(r2 >= r1, "increases", "attenuates"))

d2long <- d2 %>%
  select(gene, set, direction, `Batch only` = r1, `Fully adjusted` = r2) %>%
  pivot_longer(c(`Batch only`, `Fully adjusted`),
               names_to = "model", values_to = "ratio")

f2 <- ggplot(d2, aes(y = gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = r1, xend = r2, yend = gene, colour = direction),
               linewidth = .7,
               arrow = arrow(length = unit(0.06, "inches"), type = "closed")) +
  geom_point(data = d2long, aes(x = ratio, shape = model), size = 1.8) +
  scale_colour_manual(values = c(increases = "#1b5e20",
                                 attenuates = "#c62828"), name = NULL) +
  scale_shape_manual(values = c(`Batch only` = 1, `Fully adjusted` = 16),
                     name = NULL) +
  labs(x = "SD ratio (T2D / non-diabetic)", y = NULL,
       title = "Effect of adjusting for tissue composition",
       subtitle = "Arrows run from batch-only to fully adjusted model",
       caption = paste("Genes whose variance attenuates on adjustment had composition-sensitive signal.",
                       "\nThe crinophagy signals did not attenuate. This is consistent with independence from",
                       "\nthe measured covariates; it does not exclude unmeasured composition effects.")) +
  theme_pub

ggsave(file.path(FIGS, "Fig2_composition_adjustment.pdf"), f2,
       width = 6.5, height = 5)


# ---------------------------------------------------------------------
# FIG 3 - mean vs variance
# CROSS-COHORT. The x and y coordinates of every point come from
# different samples. This must be stated in the caption; it is the most
# likely thing a referee will challenge.
# ---------------------------------------------------------------------
# Requires the Lawlor beta-cell sheet. Skipped with a message if absent.

lawlor_file <- "Supplemental_Table_S6.xlsx"

if (file.exists(lawlor_file)) {
  library(readxl)
  beta <- read_excel(lawlor_file, sheet = "T2D.vs.NonT2D.Beta")

  d3 <- m2 %>%
    inner_join(beta %>%
                 select(gene = Associated.Gene.Name,
                        logFC, FDR_mean = FDR),
               by = "gene") %>%
    mutate(label = gene %in% c("UNC13D","STX7","CTSD","INS",
                               "STX1A","GABARAP","RAB3A"))

  f3 <- ggplot(d3, aes(logFC, sd_ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    geom_hline(yintercept = quantile(bg2, .99), linetype = "dotted",
               colour = "grey30") +
    geom_point(aes(colour = set, size = label, alpha = label)) +
    ggrepel::geom_text_repel(data = filter(d3, label), aes(label = gene),
                             size = 2.8, min.segment.length = 0,
                             box.padding = .4, seed = 1) +
    scale_colour_manual(values = COL, name = NULL) +
    scale_size_manual(values = c(`FALSE` = 1, `TRUE` = 2.4), guide = "none") +
    scale_alpha_manual(values = c(`FALSE` = .35, `TRUE` = 1), guide = "none") +
    labs(x = "Mean change in T2D beta cells (log fold change)",
         y = "Inter-donor SD ratio (T2D / non-diabetic)",
         title = "Mean and variance capture different genes",
         subtitle = "Dotted line: 99th percentile of the genome-wide variance distribution",
         caption = paste("CROSS-COHORT PLOT. The x coordinate derives from single-cell beta-cell data",
                         "\n(Lawlor et al.); the y coordinate from bulk LCM islet RNA-seq (GSE164416).",
                         "\nNo gene's two coordinates come from the same samples. Only genes tested in",
                         "\nboth datasets are shown.")) +
    theme_pub

  ggsave(file.path(FIGS, "Fig3_mean_vs_variance.pdf"), f3,
         width = 6.5, height = 5)
} else {
  message("Fig 3 skipped: ", lawlor_file, " not found in ", getwd())
}


# ---------------------------------------------------------------------
# FIG 4 - forest plot, all 46 genes
# The complete pre-registered panel, nothing hidden.
# ---------------------------------------------------------------------

f4 <- m2 %>%
  mutate(gene = reorder(gene, sd_ratio)) %>%
  ggplot(aes(sd_ratio, gene, colour = set)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = quantile(bg2, .99), linetype = "dotted",
             colour = "grey20") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), size = .25,
                  linewidth = .4) +
  scale_colour_manual(values = COL, name = NULL) +
  labs(x = "SD ratio (T2D / non-diabetic), 95% bootstrap CI", y = NULL,
       title = "Complete pre-registered panel, fully adjusted",
       subtitle = "31 crinophagy genes, 14 comparator genes, 1 post-hoc comparator",
       caption = paste("Gene sets were frozen with a timestamped DOI before analysis.",
                       "STX1A was added post-hoc\nbecause Lawlor et al. identified it as mean-shifted, and is labelled as such.")) +
  theme_pub + theme(axis.text.y = element_text(size = 7))

ggsave(file.path(FIGS, "Fig4_forest_all_genes.pdf"), f4,
       width = 6.5, height = 8)


# ---------------------------------------------------------------------
# FIG 5 - cross-platform replication
# Shows the non-replication alongside the replication. A reader who
# sees a failure reported plainly trusts the success more.
# ---------------------------------------------------------------------

rep_dat <- data.frame(
  gene     = c("UNC13D", "UNC13D", "RAB3A", "RAB3A"),
  platform = factor(c("Discovery\n(proteomics, n=123)",
                      "Replication\n(RNA-seq, n=54)",
                      "Discovery\n(proteomics, n=123)",
                      "Replication\n(RNA-seq, n=54)"),
                    levels = c("Discovery\n(proteomics, n=123)",
                               "Replication\n(RNA-seq, n=54)")),
  ratio    = c(2.065, 2.001, 1.725, 1.101),
  lo       = c(1.201, NA, 1.088, NA),
  hi       = c(2.745, NA, 2.276, NA)
)
# fill replication CIs from the frozen table
rep_dat$lo[c(2,4)] <- m2$ci_lo[match(c("UNC13D","RAB3A"), m2$gene)]
rep_dat$hi[c(2,4)] <- m2$ci_hi[match(c("UNC13D","RAB3A"), m2$gene)]

f5 <- ggplot(rep_dat, aes(platform, ratio, group = gene, colour = gene)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(linewidth = .7) +
  geom_pointrange(aes(ymin = lo, ymax = hi), size = .45) +
  scale_colour_manual(values = c(UNC13D = "#1b5e20", RAB3A = "#c62828"),
                      name = NULL) +
  labs(x = NULL, y = "SD ratio (T2D / non-diabetic), 95% CI",
       title = "Cross-platform replication",
       subtitle = "UNC13D replicates; RAB3A does not",
       caption = paste("Discovery values are per-protein variance ratios from islet proteomics;",
                       "replication values\nare from bulk RNA-seq under the fully adjusted model.",
                       "Reporting the non-replication\nalongside the replication reflects that candidates were fixed before testing.")) +
  theme_pub

ggsave(file.path(FIGS, "Fig5_cross_platform.pdf"), f5,
       width = 5.5, height = 4.5)


cat("\nWrote to", FIGS, ":\n")
print(list.files(FIGS, pattern = "\\.pdf$"))
cat("\nAll values read from frozen outputs. No new statistics computed.\n")
