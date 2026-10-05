# =====================================================================
# Supplementary figure: Rab regulatory layer
#
# WHAT THIS SHOWS. Eight Rab GEF/GAP genes examined post-hoc, plus
# RAB3GAP1 from the frozen set, against the fusion-execution genes that
# carry the signal. The regulatory layer is flat; the execution layer
# is not.
#
# WHAT IT DOES NOT SHOW. Whether the Rab27A:Rab3A activity ratio shifts
# in T2D. GEFs and GAPs act catalytically, so a change in their activity
# or localisation can alter the nucleotide state of a much larger Rab
# pool without any change in their own abundance. A flat abundance
# result is therefore uninformative about activity, not evidence
# against it. The caption says so explicitly.
#
# STATUS. Post-hoc. MADD, TBC1D10A/B, RAB3GAP2, RAB3IP, RAB3IL1,
# DENND2B and TBC1D3 were NOT in the frozen sets and were examined after
# seeing the primary results. Only RAB3GAP1, RAB3A and RAB27A were
# pre-registered. The figure distinguishes them by shape.
#
# Requires: ratio_all, bg2, mat, pd2 in session.
# =====================================================================

library(ggplot2)
library(dplyr)
library(car)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

stopifnot(exists("ratio_all"), exists("bg2"), exists("mat"), exists("pd2"))

# ---- Genes, grouped by role in the pathway --------------------------

genes <- tribble(
  ~gene,       ~ensembl,           ~layer,         ~frozen,
  # execution: these carry the signal
  "UNC13D",    "ENSG00000092929",  "Fusion execution",  TRUE,
  "STX7",      "ENSG00000079950",  "Fusion execution",  TRUE,
  "CTSD",      "ENSG00000117984",  "Lysosomal protease", TRUE,
  # Rab substrates
  "RAB27A",    "ENSG00000069974",  "Rab GTPase",        TRUE,
  "RAB3A",     "ENSG00000105649",  "Rab GTPase",        TRUE,
  # regulatory layer
  "RAB3GAP1",  "ENSG00000115839",  "Rab GEF / GAP",     TRUE,
  "MADD",      "ENSG00000110514",  "Rab GEF / GAP",     FALSE,
  "RAB3GAP2",  "ENSG00000118873",  "Rab GEF / GAP",     FALSE,
  "TBC1D10A",  NA,                 "Rab GEF / GAP",     FALSE,
  "TBC1D10B",  NA,                 "Rab GEF / GAP",     FALSE,
  "RAB3IP",    NA,                 "Rab GEF / GAP",     FALSE,
  "RAB3IL1",   NA,                 "Rab GEF / GAP",     FALSE,
  "DENND2B",   NA,                 "Rab GEF / GAP",     FALSE,
  "TBC1D3",    NA,                 "Rab GEF / GAP",     FALSE
)

# Resolve any missing accessions programmatically rather than by hand
if (any(is.na(genes$ensembl))) {
  library(org.Hs.eg.db); library(AnnotationDbi)
  need <- genes$gene[is.na(genes$ensembl)]
  mp <- AnnotationDbi::select(org.Hs.eg.db, keys = need,
                              keytype = "SYMBOL", columns = "ENSEMBL")
  mp <- mp[!is.na(mp$ENSEMBL) & mp$ENSEMBL %in% rownames(mat), ]
  mp <- mp[!duplicated(mp$SYMBOL), ]
  genes$ensembl[match(mp$SYMBOL, genes$gene)] <- mp$ENSEMBL
}

missing <- genes$gene[is.na(genes$ensembl) |
                      !genes$ensembl %in% rownames(mat)]
if (length(missing))
  cat("Not measured in this matrix:", paste(missing, collapse = ", "), "\n")

genes <- genes %>% filter(!is.na(ensembl), ensembl %in% rownames(mat))


# ---- Compute ratio, CI and percentile under the full model ----------

set.seed(1)

res <- lapply(seq_len(nrow(genes)), function(i) {
  d <- pd2
  d$y <- mat[genes$ensembl[i], rownames(d)]
  d$r <- residuals(lm(y ~ batch + imm + beta, data = d))
  a <- d$r[d$group == "T2D"]; b <- d$r[d$group == "ND"]

  boot <- replicate(5000,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))

  data.frame(gene = genes$gene[i],
             layer = genes$layer[i],
             frozen = genes$frozen[i],
             sd_ratio = sd(a) / sd(b),
             ci_lo = quantile(boot, .025),
             ci_hi = quantile(boot, .975),
             p = leveneTest(r ~ group, data = d)$`Pr(>F)`[1],
             pct = 100 * mean(bg2 <= (sd(a) / sd(b))),
             row.names = NULL)
}) %>% bind_rows()

res <- res %>%
  mutate(layer = factor(layer, levels = c("Fusion execution",
                                          "Lysosomal protease",
                                          "Rab GTPase",
                                          "Rab GEF / GAP")),
         status = ifelse(frozen, "Pre-registered", "Post-hoc")) %>%
  arrange(layer, desc(sd_ratio))

cat("\n===== RAB REGULATORY LAYER =====\n")
print(res %>% dplyr::select(gene, layer, status, sd_ratio, pct, p),
      digits = 3, row.names = FALSE)

cat("\nMedian ratio by layer:\n")
print(res %>% group_by(layer) %>%
        summarise(n = n(), median = median(sd_ratio)) %>%
        as.data.frame(), digits = 3, row.names = FALSE)

write.csv(res, "figures/rab_regulatory_layer.csv", row.names = FALSE)


# ---- Figure ---------------------------------------------------------

p99 <- quantile(bg2, .99)
p50 <- median(bg2)

ord <- res %>% arrange(layer, sd_ratio) %>% pull(gene)

fig <- res %>%
  mutate(gene = factor(gene, levels = ord)) %>%
  ggplot(aes(sd_ratio, gene)) +
  annotate("rect", xmin = quantile(bg2, .25), xmax = quantile(bg2, .75),
           ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_vline(xintercept = p50, colour = "grey55", linewidth = .3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = p99, linetype = "dotted", colour = "grey20") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi,
                      colour = layer, shape = status),
                  size = .35, linewidth = .45) +
  scale_colour_manual(values = c("Fusion execution"   = "#1b5e20",
                                 "Lysosomal protease" = "#2e7d32",
                                 "Rab GTPase"         = "#8d6e63",
                                 "Rab GEF / GAP"      = "#5c6bc0"),
                      name = NULL) +
  scale_shape_manual(values = c(`Pre-registered` = 16, `Post-hoc` = 1),
                     name = NULL) +
  facet_grid(layer ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  labs(x = "SD ratio (T2D / non-diabetic), 95% bootstrap CI", y = NULL,
       title = "Variance is concentrated in fusion execution, not Rab regulation",
       subtitle = paste0("Fully adjusted model. Grey band: interquartile range of the ",
                         format(length(bg2), big.mark = ","),
                         "-gene background.\nDashed line: no difference. Dotted line: 99th percentile."),
       caption = paste(
         "Open circles mark genes examined post-hoc; these were not in the frozen sets and are",
         "\nexcluded from all set-level statistics. GEF and GAP proteins act catalytically, so their",
         "\nabundance is a poor proxy for Rab nucleotide state: a flat result here neither supports",
         "\nnor excludes a shift in the Rab27A:Rab3A activity ratio, which would require",
         "\nnucleotide-state or effector-pulldown assays across donors.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text.y.left = element_text(angle = 0, size = 8,
                                         colour = "grey25"),
        strip.placement = "outside",
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey30"),
        plot.caption = element_text(size = 6.8, colour = "grey40",
                                    hjust = 0),
        legend.position = "top", legend.box = "horizontal")

ggsave(file.path(FIGS, "FigS_rab_regulatory_layer.pdf"), fig,
       width = 7, height = 6)

cat("\nWrote FigS_rab_regulatory_layer.pdf and rab_regulatory_layer.csv\n")
