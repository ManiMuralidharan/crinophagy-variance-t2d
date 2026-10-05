# =====================================================================
# The five missing figure panels, all as PDF
#
#   M1  Set-level permutation nulls        -> Figure 2C,D
#   M2  99th-percentile tail enrichment    -> Figure 2E
#   M3  Munc13 paralogues, quantitative    -> Figure 3B
#   M4  Composition covariates by group    -> Figure 4A
#   M5  Single-cell mean differences       -> Figure 5D
#   M6  CAPS family added to the exocytosis panel (replaces p9)
#
# All output is vector PDF with live text, editable in Illustrator.
#
# REQUIRES in session: mat, pd2, bg1, bg2
#   load("session_gse164416.RData")
# and on disk: replication_GSE164416/canonical_full.csv
#
# Run:  source("29_missing_panels.R")
# =====================================================================

library(ggplot2)
library(dplyr)
library(car)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

if (!exists("mat")) load("session_gse164416.RData")
stopifnot(exists("mat"), exists("pd2"), exists("bg1"), exists("bg2"))

set.seed(1)
NPERM <- 10000

# house style, matched to the panels you already have
COL  <- c(crinophagy = "#1b5e20", comparator = "#90a4ae",
          `post-hoc` = "#ef6c00")
COLG <- c(`Non-diabetic` = "#90a4ae", `Pre-T2D` = "#66bb6a",
          `T2D` = "#1b5e20")

th <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        plot.caption = element_text(size = 6.8, colour = "grey40",
                                    hjust = 0))

full <- read.csv(file.path("replication_GSE164416", "canonical_full.csv"),
                 stringsAsFactors = FALSE)

cat("Session OK. Background sizes:", length(bg1), length(bg2), "\n\n")


# =====================================================================
# M1  Set-level permutation nulls  (Figure 2C,D)
# =====================================================================

cat("=== M1: permutation nulls ===\n")

perm_null <- function(tag, bg) {
  res <- full %>% dplyr::filter(model == tag)
  out <- lapply(c("crinophagy", "comparator"), function(s) {
    g <- res$sd_ratio[res$set == s]; g <- g[is.finite(g)]
    obs <- median(g)
    null <- replicate(NPERM, median(sample(bg, length(g))))
    p <- (sum(null >= obs) + 1) / (NPERM + 1)
    list(n = data.frame(set = s, null = null),
         o = data.frame(set = s, obs = obs, n_genes = length(g), p = p))
  })
  list(nulls = bind_rows(lapply(out, `[[`, "n")),
       obs   = bind_rows(lapply(out, `[[`, "o")))
}

m2 <- perm_null("M2_full", bg2)
print(m2$obs, digits = 4, row.names = FALSE)

lab <- c(crinophagy = "Crinophagy set (31 genes, pre-registered)",
         comparator = "Macroautophagy comparator (14 genes, pre-registered)")
m2$nulls$set <- factor(lab[m2$nulls$set], levels = lab)
m2$obs$set   <- factor(lab[m2$obs$set],   levels = lab)
m2$obs$txt <- sprintf("observed %.3f\npermutation P = %.3f",
                      m2$obs$obs, m2$obs$p)

p1 <- ggplot(m2$nulls, aes(null)) +
  geom_histogram(bins = 60, fill = "grey86", colour = NA) +
  geom_vline(data = m2$obs, aes(xintercept = obs),
             colour = "#1b5e20", linewidth = .8) +
  geom_text(data = m2$obs, aes(x = obs, y = Inf, label = txt),
            hjust = -0.08, vjust = 1.4, size = 2.7, colour = "#1b5e20") +
  facet_wrap(~ set, nrow = 1, scales = "free") +
  labs(x = "Median SD ratio of a random gene set of the same size",
       y = "Permutations",
       title = "Set-level permutation test",
       subtitle = sprintf("%s random size-matched gene sets drawn from the genome-wide reference distribution; fully adjusted model",
                          format(NPERM, big.mark = ",")),
       caption = "Both sets were frozen with a timestamped DOI before analysis, so these tests are pre-specified.") +
  th + theme(strip.text = element_text(size = 8.5, face = "bold"))

ggsave(file.path(FIGS, "M1_permutation_nulls.pdf"), p1,
       width = 8, height = 3.6, useDingbats = FALSE)
cat("  wrote M1_permutation_nulls.pdf\n")


# =====================================================================
# M2  Tail enrichment  (Figure 2E)
# =====================================================================

cat("\n=== M2: tail enrichment ===\n")

p99 <- quantile(bg2, .99)
r2  <- full %>% dplyr::filter(model == "M2_full")

tail_df <- bind_rows(lapply(c("crinophagy","comparator"), function(s) {
  g <- r2$sd_ratio[r2$set == s]; g <- g[is.finite(g)]
  k <- sum(g >= p99); n <- length(g)
  data.frame(set = s, n = n, observed = k, expected = 0.01 * n,
             p = binom.test(k, n, 0.01, alternative = "greater")$p.value,
             genes = paste(r2$gene[r2$set == s & r2$sd_ratio >= p99],
                           collapse = ", "))
}))
print(tail_df, digits = 3, row.names = FALSE)

td <- tail_df %>%
  tidyr::pivot_longer(c(observed, expected),
                      names_to = "what", values_to = "count") %>%
  mutate(set = factor(lab[set], levels = lab),
         what = factor(what, levels = c("expected","observed"),
                       labels = c("Expected under the reference distribution",
                                  "Observed")))

ann <- tail_df %>%
  mutate(set = factor(lab[set], levels = lab),
         txt = sprintf("binomial P = %.3f%s", p,
                       ifelse(nchar(genes) > 0,
                              paste0("\n", genes), "")))

p2 <- ggplot(td, aes(set, count, fill = what)) +
  geom_col(position = position_dodge(width = .65), width = .6) +
  geom_text(data = ann, aes(x = set, y = 2.35, label = txt),
            inherit.aes = FALSE, size = 2.7, colour = "grey25") +
  scale_fill_manual(values = c("grey80", "#1b5e20"), name = NULL) +
  ylim(0, 2.8) +
  labs(x = NULL, y = "Genes above the 99th percentile",
       title = "Tail enrichment at the 99th percentile",
       subtitle = sprintf("Fully adjusted model; threshold = %.3f", p99),
       caption = "Counts are of pre-registered set members exceeding the model-specific 99th percentile of the 20,255-gene reference distribution.") +
  th + theme(legend.position = "top",
             axis.text.x = element_text(size = 8))

ggsave(file.path(FIGS, "M2_tail_enrichment.pdf"), p2,
       width = 7, height = 4, useDingbats = FALSE)
cat("  wrote M2_tail_enrichment.pdf\n")


# =====================================================================
# M3  Munc13 paralogues, quantitative  (Figure 3B)
# =====================================================================

cat("\n=== M3: Munc13 paralogues ===\n")

one <- function(ens, nboot = 5000) {
  d <- pd2
  d$y <- mat[ens, rownames(d)]
  d$r <- residuals(lm(y ~ batch + imm + beta, data = d))
  a <- d$r[d$group == "T2D"]; b <- d$r[d$group == "ND"]
  bt <- replicate(nboot,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))
  list(stat = data.frame(ratio = sd(a)/sd(b),
                         lo = unname(quantile(bt, .025)),
                         hi = unname(quantile(bt, .975)),
                         p  = leveneTest(r ~ group, data = d)$`Pr(>F)`[1],
                         pct = 100 * mean(bg2 <= sd(a)/sd(b))),
       donors = data.frame(group = d$group, resid = d$r))
}

M13 <- c(`Munc13-1\n(UNC13A)` = "ENSG00000130477",
         `Munc13-2\n(UNC13B)` = "ENSG00000198722",
         `Munc13-4\n(UNC13D)` = "ENSG00000092929")
M13 <- M13[M13 %in% rownames(mat)]
if (length(M13) < 3) cat("  WARNING: only", length(M13), "paralogues found\n")

o <- Map(one, M13)
S <- bind_rows(lapply(names(o), function(n)
  cbind(label = n, o[[n]]$stat)))
D <- bind_rows(lapply(names(o), function(n)
  cbind(label = n, o[[n]]$donors)))

S$label <- factor(S$label, levels = names(M13))
D$label <- factor(D$label, levels = names(M13))
D$group <- factor(D$group, levels = c("ND","T2D"),
                  labels = c("Non-diabetic","T2D"))

cat("\n"); print(S, digits = 4, row.names = FALSE)

pa <- ggplot(D, aes(group, resid, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .55, width = .55,
               linewidth = .3) +
  geom_jitter(width = .15, size = .9, alpha = .6, colour = "grey20") +
  facet_wrap(~ label, nrow = 1) +
  scale_fill_manual(values = c(`Non-diabetic` = "#90a4ae",
                               `T2D` = "#1b5e20"), guide = "none") +
  labs(x = NULL, y = "Adjusted expression (residual)",
       title = "A   Donor-level distributions") +
  th + theme(strip.text = element_text(size = 8, lineheight = .95))

pb <- ggplot(S, aes(ratio, label, colour = label)) +
  annotate("rect", xmin = quantile(bg2, .25), xmax = quantile(bg2, .75),
           ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_vline(xintercept = median(bg2), colour = "grey55", linewidth = .3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = p99, linetype = "dotted", colour = "grey20") +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = .45, linewidth = .6) +
  geom_text(aes(x = hi, label = sprintf("  %.1fth pct", pct)),
            hjust = 0, size = 2.6, colour = "grey30") +
  scale_colour_manual(values = c("#5c6bc0","#8d6e63","#1b5e20"),
                      guide = "none") +
  xlim(0.3, 3.6) +
  labs(x = "SD ratio (T2D / non-diabetic), 95% bootstrap CI", y = NULL,
       title = "B   Effect sizes against the genome-wide background",
       caption = paste("Munc13-1 and Munc13-2 prime secretory granules for fusion with the plasma membrane;",
                       "\nMunc13-4 primes them for fusion with lysosomes. Grey band: interquartile range of the",
                       "\n20,255-gene background. Munc13-3 (UNC13C) was not detected in islet tissue.")) +
  th + theme(panel.grid.major.y = element_blank(),
             axis.text.y = element_text(size = 8, lineheight = .95))

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  ggsave(file.path(FIGS, "M3_munc13_paralogues.pdf"),
         pa / pb + plot_layout(heights = c(1, 0.9)),
         width = 6.5, height = 6, useDingbats = FALSE)
} else {
  ggsave(file.path(FIGS, "M3_munc13_A.pdf"), pa, width = 6.5, height = 3)
  ggsave(file.path(FIGS, "M3_munc13_B.pdf"), pb, width = 6.5, height = 3)
}
cat("  wrote M3_munc13_paralogues.pdf\n")


# =====================================================================
# M4  Composition covariates by group  (Figure 4A)
# =====================================================================

cat("\n=== M4: composition covariates ===\n")

cv <- bind_rows(
  data.frame(group = pd2$group, value = pd2$imm,
             score = "Immune infiltration score"),
  data.frame(group = pd2$group, value = pd2$beta,
             score = "Beta-cell content score")) %>%
  mutate(group = factor(group, levels = c("ND","T2D"),
                        labels = c("Non-diabetic","T2D")),
         score = factor(score, levels = c("Immune infiltration score",
                                          "Beta-cell content score")))

sds <- cv %>% group_by(score, group) %>%
  summarise(sd = sd(value), n = n(), .groups = "drop")
print(as.data.frame(sds), digits = 3, row.names = FALSE)

ann4 <- sds %>% group_by(score) %>%
  summarise(txt = sprintf("SD %.3f vs %.3f",
                          sd[group == "T2D"], sd[group == "Non-diabetic"]),
            .groups = "drop")

p4 <- ggplot(cv, aes(group, value, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .55, width = .55,
               linewidth = .3) +
  geom_jitter(width = .15, size = 1, alpha = .6, colour = "grey20") +
  geom_text(data = ann4, aes(x = 1.5, y = Inf, label = txt),
            inherit.aes = FALSE, vjust = 1.6, size = 2.8,
            colour = "grey25") +
  facet_wrap(~ score, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = c(`Non-diabetic` = "#90a4ae",
                               `T2D` = "#1b5e20"), guide = "none") +
  labs(x = NULL, y = "Score (mean variance-stabilised expression)",
       title = "Tissue composition varies more among T2D donors",
       subtitle = "Replication cohort, 18 non-diabetic and 36 T2D donors",
       caption = paste("Both composition scores are more variable among T2D donors, which is why adjustment",
                       "\nwas necessary rather than optional. Immune score: 34 marker genes. Beta-cell score: 8.")) +
  th + theme(strip.text = element_text(size = 8.5, face = "bold"))

ggsave(file.path(FIGS, "M4_composition_covariates.pdf"), p4,
       width = 6.5, height = 4, useDingbats = FALSE)
cat("  wrote M4_composition_covariates.pdf\n")


# =====================================================================
# M5  Single-cell mean differences  (Figure 5D)
# =====================================================================

cat("\n=== M5: single-cell means ===\n")

# Values from Lawlor et al. supplementary table S6, sheet
# T2D.vs.NonT2D.Beta. Hard-coded because no reanalysis was performed;
# these are published values, not computed here.
sc <- data.frame(
  gene  = c("INS","UNC13D","STX7","CTSD"),
  logFC = c(-0.540, 0.149, -0.392, -0.296),
  fdr   = c(0.0060, 0.741,  0.894,  0.838),
  role  = c("Cargo","Machinery","Machinery","Machinery"))
sc$gene <- factor(sc$gene, levels = rev(c("INS","UNC13D","STX7","CTSD")))
sc$sig  <- ifelse(sc$fdr < 0.05, "FDR < 0.05", "not significant")
print(sc, row.names = FALSE)

p5 <- ggplot(sc, aes(logFC, gene, fill = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_col(width = .55) +
  geom_text(aes(x = logFC + ifelse(logFC < 0, -0.03, 0.03),
                label = sprintf("FDR %.3f", fdr),
                hjust = ifelse(logFC < 0, 1, 0)),
            size = 2.7, colour = "grey25") +
  scale_fill_manual(values = c(`FDR < 0.05` = "#c62828",
                               `not significant` = "grey78"),
                    name = NULL) +
  xlim(-0.85, 0.5) +
  labs(x = "Mean difference in T2D beta cells (log fold change)", y = NULL,
       title = "No mean shift in the machinery; insulin cargo is reduced",
       subtitle = "Published single-cell beta-cell data (Lawlor et al.)",
       caption = paste("Values are from the published supplementary tables; no reanalysis was performed. The",
                       "\nsame analysis detected 248 of 7,441 genes at FDR below 0.05, indicating adequate power.")) +
  th + theme(legend.position = "top",
             panel.grid.major.y = element_blank())

ggsave(file.path(FIGS, "M5_singlecell_means.pdf"), p5,
       width = 6, height = 3.6, useDingbats = FALSE)
cat("  wrote M5_singlecell_means.pdf\n")


# =====================================================================
# M6  Exocytosis panel WITH the CAPS family  (replaces p9)
# =====================================================================

cat("\n=== M6: exocytosis panel including CAPS family ===\n")

panel <- data.frame(
  gene = c("STX7","CTSD",
           "UNC13A","UNC13B","UNC13D",
           "CADPS","CADPS2","BAIAP3",
           "STX1A","STX4","SNAP25","SNAP23","VAMP2",
           "STXBP1","STXBP2","STXBP3",
           "STXBP5","STXBP5L",
           "SYT7","SYT9","RIMS2"),
  group = c(rep("Crinophagy reference", 2),
            rep("Munc13 family (MUN domain)", 3),
            rep("CAPS family (MUN domain)", 3),
            rep("Plasma membrane SNARE", 5),
            rep("Munc18 / SM protein", 3),
            rep("Tomosyn", 2),
            rep("Ca2+ sensor / active zone", 3)),
  stringsAsFactors = FALSE)

frozen <- c("UNC13D","UNC13B","STX7","CTSD","SNAP23","CADPS")

library(org.Hs.eg.db); library(AnnotationDbi)
mp <- AnnotationDbi::select(org.Hs.eg.db, keys = panel$gene,
                            keytype = "SYMBOL", columns = "ENSEMBL")
mp <- mp[!is.na(mp$ENSEMBL) & mp$ENSEMBL %in% rownames(mat), ]
mp <- mp[!duplicated(mp$SYMBOL), ]
detach("package:org.Hs.eg.db", unload = TRUE)
detach("package:AnnotationDbi", unload = TRUE)

miss <- setdiff(panel$gene, mp$SYMBOL)
if (length(miss)) cat("  not measured:", paste(miss, collapse = ", "), "\n")
panel <- merge(panel, mp, by.x = "gene", by.y = "SYMBOL")
names(panel)[names(panel) == "ENSEMBL"] <- "ensembl"

res6 <- do.call(rbind, lapply(seq_len(nrow(panel)), function(i) {
  s <- one(panel$ensembl[i])$stat
  cbind(gene = panel$gene[i], group = panel$group[i], s)
}))
res6$status <- ifelse(res6$gene %in% frozen, "Pre-registered", "Post-hoc")
res6 <- res6[order(-res6$ratio), ]

cat("\n"); print(res6[, c("gene","group","status","ratio","pct","p")],
                 digits = 3, row.names = FALSE)

exo <- res6[res6$group != "Crinophagy reference", ]
cat(sprintf("\nExocytosis genes above the 99th percentile: %d of %d\n",
            sum(exo$ratio >= p99), nrow(exo)))

glev <- c("Crinophagy reference","Munc13 family (MUN domain)",
          "CAPS family (MUN domain)","Plasma membrane SNARE",
          "Munc18 / SM protein","Tomosyn","Ca2+ sensor / active zone")
res6$group <- factor(res6$group, levels = glev)
ordg <- res6[order(res6$group, res6$ratio), "gene"]
res6$gene <- factor(res6$gene, levels = ordg)

p6 <- ggplot(res6, aes(ratio, gene)) +
  annotate("rect", xmin = quantile(bg2, .25), xmax = quantile(bg2, .75),
           ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_vline(xintercept = median(bg2), colour = "grey55", linewidth = .3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = p99, linetype = "dotted", colour = "grey20") +
  geom_pointrange(aes(xmin = lo, xmax = hi, colour = group,
                      shape = status), size = .34, linewidth = .45) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_colour_manual(values = c("#1b5e20","#ef6c00","#8d6e63",
                                 "#5c6bc0","#00838f","#7b1fa2",
                                 "#c2185b"), guide = "none") +
  scale_shape_manual(values = c(`Pre-registered` = 16, `Post-hoc` = 1),
                     name = NULL) +
  labs(x = "SD ratio (T2D / non-diabetic), 95% bootstrap CI", y = NULL,
       title = "Priming factors and exocytosis machinery as specificity controls",
       subtitle = "Fully adjusted model. Grey band: interquartile range of the 20,255-gene background. Dotted line: 99th percentile.",
       caption = paste("POST-HOC PANEL, excluded from all set-level statistics. Munc13 and CAPS proteins both carry",
                       "\nMUN domains and catalyse the same SNARE-opening step. Of the MUN-domain priming factors",
                       "\ndetected in islet tissue, only the lysosome-directed Munc13-4 shows elevated variance.",
                       "\nOpen circles: post-hoc. Filled circles: pre-registered set members shown for reference.")) +
  th + theme(panel.grid.major.y = element_blank(),
             axis.text.y = element_text(size = 7.5),
             strip.text.y.left = element_text(angle = 0, size = 7,
                                              colour = "grey25"),
             strip.placement = "outside",
             legend.position = "top")

ggsave(file.path(FIGS, "M6_exocytosis_with_CAPS.pdf"), p6,
       width = 7.5, height = 7.5, useDingbats = FALSE)
cat("  wrote M6_exocytosis_with_CAPS.pdf\n")


# ---- summary ---------------------------------------------------------

write.csv(res6, file.path(FIGS, "M6_panel_values.csv"), row.names = FALSE)
write.csv(S,    file.path(FIGS, "M3_munc13_values.csv"), row.names = FALSE)
write.csv(tail_df, file.path(FIGS, "M2_tail_values.csv"), row.names = FALSE)

cat("\n===== DONE =====\n")
print(list.files(FIGS, pattern = "^M[0-9]_.*\\.pdf$"))
cat("\nAll vector PDFs with live text. Assemble in Illustrator.\n")
