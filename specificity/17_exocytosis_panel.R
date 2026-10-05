# =====================================================================
# Post-hoc specificity panel: regulated exocytosis machinery
#
# PURPOSE, DECLARED BEFORE RUNNING. The T2D islet literature centres on
# the regulated exocytosis machinery: STX1A, SNAP25, VAMP2, Munc18
# isoforms, Munc13 isoforms, tomosyn, synaptotagmins, RIMS2. If those
# genes show ordinary inter-donor variance while the crinophagy fusion
# components are extreme, the phenotype is not a general property of
# beta-cell membrane fusion. That is the control this panel provides.
#
# The rationale comes from published biology, not from the observed
# ranks. That is what distinguishes it from a fishing expedition.
#
# INTERPRETATION FIXED IN ADVANCE:
#
#   (a) Exocytosis genes mostly ordinary (median near the background
#       median, none or one above the 99th percentile)
#       => the variance phenotype is specific to the crinophagy fusion
#          arm and not a general feature of beta-cell fusion machinery.
#
#   (b) Several exocytosis genes extreme
#       => the phenotype is BROADER than crinophagy. Report that. Do
#          not reframe the panel, drop genes, or reclassify anything to
#          preserve the crinophagy claim. UNC13B is already known to be
#          flat (0.821, 9.7th percentile), so a broad result would not
#          be uniform, and that nuance would need reporting too.
#
# STATUS. POST-HOC. Not in the frozen sets. Excluded from every
# set-level statistic. Must be declared as post-hoc in the Methods and
# marked as such in the figure.
#
# Requires: mat, pd2, bg2 in session.
# =====================================================================

library(ggplot2)
library(dplyr)
library(car)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

stopifnot(exists("mat"), exists("pd2"), exists("bg2"))
set.seed(1)


# ---- Panel definition, by function ----------------------------------

panel <- tribble(
  ~gene,      ~group,
  # plasma membrane SNAREs
  "STX1A",    "Plasma membrane SNARE",
  "STX1B",    "Plasma membrane SNARE",
  "STX4",     "Plasma membrane SNARE",
  "SNAP25",   "Plasma membrane SNARE",
  "SNAP23",   "Plasma membrane SNARE",
  "VAMP2",    "Plasma membrane SNARE",
  # SM proteins
  "STXBP1",   "Munc18 / SM protein",
  "STXBP2",   "Munc18 / SM protein",
  "STXBP3",   "Munc18 / SM protein",
  # tomosyn: inhibitory SNARE regulator
  "STXBP5",   "Tomosyn",
  "STXBP5L",  "Tomosyn",
  # Munc13 family
  "UNC13A",   "Munc13 family",
  "UNC13B",   "Munc13 family",
  "UNC13C",   "Munc13 family",
  "UNC13D",   "Munc13 family",
  # Ca2+ sensors and active zone
  "SYT7",     "Ca2+ sensor / active zone",
  "SYT9",     "Ca2+ sensor / active zone",
  "RIMS2",    "Ca2+ sensor / active zone",
  # crinophagy reference points (pre-registered, for context)
  "STX7",     "Crinophagy reference",
  "CTSD",     "Crinophagy reference"
)

frozen_genes <- c("UNC13D","UNC13B","STX7","CTSD","SNAP23","VAMP8")

mp <- AnnotationDbi::select(org.Hs.eg.db, keys = panel$gene,
                            keytype = "SYMBOL", columns = "ENSEMBL")
mp <- mp[!is.na(mp$ENSEMBL) & mp$ENSEMBL %in% rownames(mat), ]
mp <- mp[!duplicated(mp$SYMBOL), ]

panel <- panel %>%
  left_join(mp, by = c("gene" = "SYMBOL")) %>%
  rename(ensembl = ENSEMBL)

absent <- panel$gene[is.na(panel$ensembl)]
if (length(absent))
  cat("Not measured in this matrix:", paste(absent, collapse = ", "), "\n")
panel <- panel %>% dplyr::filter(!is.na(ensembl))

cat("Genes tested:", nrow(panel), "\n\n")


# ---- Same model, same background ------------------------------------

res <- lapply(seq_len(nrow(panel)), function(i) {
  d <- pd2
  d$y <- mat[panel$ensembl[i], rownames(d)]
  d$r <- residuals(lm(y ~ batch + imm + beta, data = d))
  a <- d$r[d$group == "T2D"]; b <- d$r[d$group == "ND"]

  boot <- replicate(5000,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))

  data.frame(gene = panel$gene[i], group = panel$group[i],
             sd_ratio = sd(a) / sd(b),
             ci_lo = quantile(boot, .025),
             ci_hi = quantile(boot, .975),
             p = leveneTest(r ~ group, data = d)$`Pr(>F)`[1],
             pct = 100 * mean(bg2 <= (sd(a) / sd(b))),
             row.names = NULL)
}) %>% bind_rows()

res$frozen <- res$gene %in% frozen_genes
res$status <- ifelse(res$frozen, "Pre-registered", "Post-hoc")

res <- res %>% arrange(desc(sd_ratio))

cat("===== EXOCYTOSIS SPECIFICITY PANEL =====\n")
print(res %>% dplyr::select(gene, group, status, sd_ratio, pct, p),
      digits = 3, row.names = FALSE)


# ---- Apply the pre-declared decision rule ---------------------------

p99 <- quantile(bg2, .99)
exo <- res %>% dplyr::filter(group != "Crinophagy reference")
n_hi <- sum(exo$sd_ratio >= p99)
med_exo <- median(exo$sd_ratio)
med_bg  <- median(bg2)

cat(sprintf("\n99th percentile threshold: %.3f\n", p99))
cat(sprintf("Exocytosis panel: n=%d, median=%.3f (background median %.3f)\n",
            nrow(exo), med_exo, med_bg))
cat(sprintf("Above 99th percentile: %d (expect %.2f by chance)\n",
            n_hi, 0.01 * nrow(exo)))

cat("\nVERDICT (rule fixed before running):\n")
if (n_hi <= 1 && med_exo < 1.25) {
  cat("  (a) The exocytosis machinery is largely ordinary. The variance\n")
  cat("      phenotype is not a general feature of beta-cell membrane\n")
  cat("      fusion. This supports specificity to the crinophagy arm.\n")
} else {
  cat("  (b) Multiple exocytosis genes are extreme. The phenotype is\n")
  cat("      BROADER than crinophagy and must be reported as such.\n")
  cat("      Do not reframe the panel or drop genes. Note also that\n")
  cat("      UNC13B sits at the 9.7th percentile, so any broad effect\n")
  cat("      is not uniform across the machinery.\n")
}

write.csv(res, file.path(FIGS, "exocytosis_panel.csv"), row.names = FALSE)


# ---- Figure ---------------------------------------------------------

ord <- res %>% arrange(group, sd_ratio) %>% pull(gene)

fig <- res %>%
  mutate(gene = factor(gene, levels = ord),
         group = factor(group, levels = c("Crinophagy reference",
                                          "Munc13 family",
                                          "Plasma membrane SNARE",
                                          "Munc18 / SM protein",
                                          "Tomosyn",
                                          "Ca2+ sensor / active zone"))) %>%
  ggplot(aes(sd_ratio, gene)) +
  annotate("rect", xmin = quantile(bg2, .25), xmax = quantile(bg2, .75),
           ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_vline(xintercept = median(bg2), colour = "grey55", linewidth = .3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = p99, linetype = "dotted", colour = "grey20") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi,
                      colour = group, shape = status),
                  size = .35, linewidth = .45) +
  scale_colour_manual(values = c("Crinophagy reference"      = "#1b5e20",
                                 "Munc13 family"             = "#ef6c00",
                                 "Plasma membrane SNARE"     = "#5c6bc0",
                                 "Munc18 / SM protein"       = "#00838f",
                                 "Tomosyn"                   = "#8d6e63",
                                 "Ca2+ sensor / active zone" = "#7b1fa2"),
                      name = NULL) +
  scale_shape_manual(values = c(`Pre-registered` = 16, `Post-hoc` = 1),
                     name = NULL) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  guides(colour = "none") +
  labs(x = "SD ratio (T2D / non-diabetic), 95% bootstrap CI", y = NULL,
       title = "Regulated exocytosis machinery as a specificity control",
       subtitle = paste0("Fully adjusted model. Grey band: interquartile range of the ",
                         format(length(bg2), big.mark = ","),
                         "-gene background.\nDotted line: 99th percentile."),
       caption = paste(
         "POST-HOC PANEL. These genes were not in the pre-registered sets and are excluded from all",
         "\nset-level statistics. They were selected from the published T2D insulin-secretion literature",
         "\nto test whether elevated inter-donor variance is a general property of beta-cell membrane",
         "\nfusion machinery. Open circles: post-hoc. Filled circles: pre-registered set members shown",
         "\nfor reference.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 7.5,
                                         colour = "grey25"),
        strip.placement = "outside",
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        plot.caption = element_text(size = 6.8, colour = "grey40", hjust = 0),
        legend.position = "top")

ggsave(file.path(FIGS, "FigS_exocytosis_specificity.pdf"), fig,
       width = 7, height = 7)

cat("\nWrote FigS_exocytosis_specificity.pdf and exocytosis_panel.csv\n")
