# =====================================================================
# Munc13 paralogue figure
#
# THE POINT. Three paralogues of one gene family, defined by sequence
# homology rather than by any grouping we chose after seeing results:
#
#   UNC13D (Munc13-4)  2.001  99.4th percentile
#   UNC13A (Munc13-1)  1.086  51.5th
#   UNC13B (Munc13-2)  0.763   5.9th
#
# Munc13-1 and Munc13-2 prime secretory granules for fusion with the
# plasma membrane (exocytosis). Munc13-4 primes them for fusion with
# lysosomes (crinophagy). The paralogue that diverges in T2D is the
# degradation-directed one.
#
# This is the cleanest specificity argument in the paper because the
# grouping is a gene family, not a functional partition we drew.
#
# UNC13C (Munc13-3) is largely cerebellar and was not detected.
#
# Requires: mat, pd2, bg2 in session.
# =====================================================================

library(ggplot2)
library(dplyr)
library(car)
library(patchwork)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

stopifnot(exists("mat"), exists("pd2"), exists("bg2"))
set.seed(1)

M13 <- c(`Munc13-1\n(UNC13A)` = "ENSG00000130477",
         `Munc13-2\n(UNC13B)` = "ENSG00000198722",
         `Munc13-4\n(UNC13D)` = "ENSG00000092929")

# Resolve programmatically rather than trusting the IDs above
if (!all(M13 %in% rownames(mat))) {
  library(org.Hs.eg.db); library(AnnotationDbi)
  mp <- AnnotationDbi::select(org.Hs.eg.db,
                              keys = c("UNC13A","UNC13B","UNC13D"),
                              keytype = "SYMBOL", columns = "ENSEMBL")
  mp <- mp[!is.na(mp$ENSEMBL) & mp$ENSEMBL %in% rownames(mat), ]
  mp <- mp[!duplicated(mp$SYMBOL), ]
  M13 <- setNames(mp$ENSEMBL,
                  c(`UNC13A` = "Munc13-1\n(UNC13A)",
                    `UNC13B` = "Munc13-2\n(UNC13B)",
                    `UNC13D` = "Munc13-4\n(UNC13D)")[mp$SYMBOL])
  detach("package:org.Hs.eg.db", unload = TRUE)
  detach("package:AnnotationDbi", unload = TRUE)
}

cat("Paralogues found:", length(M13), "\n")


# ---- Stats ----------------------------------------------------------

stat1 <- function(ens, label) {
  d <- pd2
  d$y <- mat[ens, rownames(d)]
  d$r <- residuals(lm(y ~ batch + imm + beta, data = d))
  a <- d$r[d$group == "T2D"]; b <- d$r[d$group == "ND"]
  boot <- replicate(5000,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))
  list(
    summary = data.frame(
      label = label, ensembl = ens,
      sd_ratio = sd(a) / sd(b),
      ci_lo = quantile(boot, .025), ci_hi = quantile(boot, .975),
      p = leveneTest(r ~ group, data = d)$`Pr(>F)`[1],
      pct = 100 * mean(bg2 <= (sd(a) / sd(b))),
      sd_T2D = sd(a), sd_ND = sd(b), row.names = NULL),
    donors = data.frame(label = label, group = d$group, resid = d$r)
  )
}

out <- Map(stat1, M13, names(M13))
S <- bind_rows(lapply(out, `[[`, "summary"))
D <- bind_rows(lapply(out, `[[`, "donors"))

S$label <- factor(S$label, levels = names(M13))
D$label <- factor(D$label, levels = names(M13))
D$group <- factor(D$group, levels = c("ND", "T2D"),
                  labels = c("Non-diabetic", "T2D"))

cat("\n===== MUNC13 PARALOGUES =====\n")
print(S %>% dplyr::select(label, sd_ratio, ci_lo, ci_hi, p, pct,
                          sd_ND, sd_T2D),
      digits = 3, row.names = FALSE)

write.csv(S, file.path(FIGS, "munc13_paralogues.csv"), row.names = FALSE)

p99 <- quantile(bg2, .99)
p50 <- median(bg2)


# ---- Panel A: donor-level distributions -----------------------------
# Shows what a variance difference actually looks like: same centre,
# different spread. Readers used to mean comparisons need to see this.

pa <- ggplot(D, aes(group, resid, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .55, width = .55,
               linewidth = .3) +
  geom_jitter(width = .16, size = 1, alpha = .65, colour = "grey20") +
  facet_wrap(~ label, nrow = 1) +
  scale_fill_manual(values = c(`Non-diabetic` = "#90a4ae",
                               `T2D` = "#1b5e20"), guide = "none") +
  labs(x = NULL, y = "Adjusted expression (residual)",
       title = "A   Donor-level distributions") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8.5, lineheight = .95),
        plot.title = element_text(size = 10, face = "bold"))


# ---- Panel B: ratios against the genome-wide background -------------

pb <- ggplot() +
  geom_histogram(data = data.frame(r = as.numeric(bg2)), aes(r),
                 bins = 110, fill = "grey88", colour = NA) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = p99, linetype = "dotted", colour = "grey20") +
  geom_vline(data = S, aes(xintercept = sd_ratio, colour = label),
             linewidth = .7) +
  geom_text(data = S,
            aes(x = sd_ratio, y = Inf,
                label = sub("\n", " ", label), colour = label),
            angle = 90, hjust = 1.05, vjust = -0.4, size = 2.5) +
  scale_colour_manual(values = c("#5c6bc0", "#8d6e63", "#1b5e20"),
                      guide = "none") +
  coord_cartesian(xlim = c(0.3, 2.6)) +
  labs(x = "SD ratio (T2D / non-diabetic)", y = "Genes",
       title = "B   Position in the genome-wide distribution",
       subtitle = sprintf("%s genes. Dotted line: 99th percentile (%.2f)",
                          format(length(bg2), big.mark = ","), p99)) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 7.5, colour = "grey35"))


# ---- Panel C: ratios with CIs ---------------------------------------

pc <- ggplot(S, aes(sd_ratio, label, colour = label)) +
  annotate("rect", xmin = quantile(bg2, .25), xmax = quantile(bg2, .75),
           ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_vline(xintercept = p50, colour = "grey55", linewidth = .3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = p99, linetype = "dotted", colour = "grey20") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), size = .45,
                  linewidth = .6) +
  geom_text(aes(x = ci_hi, label = sprintf("  %.1fth pct", pct)),
            hjust = 0, size = 2.6, colour = "grey30") +
  scale_colour_manual(values = c("#5c6bc0", "#8d6e63", "#1b5e20"),
                      guide = "none") +
  xlim(0.3, 3.6) +
  labs(x = "SD ratio (T2D / non-diabetic), 95% bootstrap CI", y = NULL,
       title = "C   Effect sizes",
       caption = paste(
         "Munc13-1 and Munc13-2 prime secretory granules for fusion with the plasma membrane;",
         "\nMunc13-4 primes them for fusion with lysosomes. Grey band: interquartile range of the",
         "\ngenome-wide background. Munc13-3 (UNC13C) is largely cerebellar and was not detected.",
         "\nModel adjusted for batch, immune infiltration and beta-cell content.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 8, lineheight = .95),
        plot.title = element_text(size = 10, face = "bold"),
        plot.caption = element_text(size = 6.8, colour = "grey40", hjust = 0))


fig <- pa / pb / pc + plot_layout(heights = c(1.15, 1, 0.85))

ggsave(file.path(FIGS, "Fig_munc13_paralogues.pdf"), fig,
       width = 6.5, height = 9)

cat("\nWrote Fig_munc13_paralogues.pdf and munc13_paralogues.csv\n")
cat("\nIf patchwork is unavailable: install.packages('patchwork'),\n")
cat("or save pa, pb, pc separately and assemble in Illustrator.\n")
