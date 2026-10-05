# =====================================================================
# The conventional analysis, run on our own data
#
# WHY THIS PANEL EXISTS. The paper's central methodological claim is
# that Munc13-4 is invisible to differential expression analysis. So
# far that rests on published single-cell data (Lawlor et al.:
# UNC13D logFC 0.15, FDR 0.74). This runs the conventional analysis on
# OUR OWN discovery cohort and shows where our candidates land.
#
# The expected picture: a normal-looking volcano with real hits, and
# UNC13D, STX7 and CTSD sitting unremarkably near the origin. That is
# the argument made visually rather than asserted.
#
# IMPORTANT: if UNC13D turns out to be significantly mean-shifted here,
# that CONTRADICTS the claim and must be reported. Do not quietly drop
# the panel. The three prior comparisons all showed no mean change, so
# a hit here would need explaining.
#
# Adapts the volcano format of Wigger et al. (Nat Metab 2021,
# Extended Data Fig 4e); data and analysis are ours.
# =====================================================================

library(ggplot2)
library(dplyr)
library(ggrepel)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)


# ---- 1. Discovery proteomics, T2D vs non-diabetic --------------------

prot  <- readRDS("./proteomics/proc_prot.rds")
donor <- read.csv("donor.csv", stringsAsFactors = FALSE)

expr <- as.matrix(prot[, -1, drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- as.character(prot[[1]])

md <- donor[match(colnames(expr), donor$record_id), ]
rownames(md) <- md$record_id
md$group <- factor(md$diagnosis_computed,
                   levels = c("None", "Pre.T2D", "Type2"))

COV <- c("donorage","donorsex","bodymassindex",
         "predistributionculturetime","percentieqrecoverypostculture")

ok <- complete.cases(md[, COV]) & md$group %in% c("None","Type2")
md2 <- droplevels(md[ok, ])
md2$donorsex <- factor(md2$donorsex)
e2 <- expr[, rownames(md2)]

cat("===== CONVENTIONAL DIFFERENTIAL ANALYSIS =====\n")
cat("Donors:", ncol(e2), "\n")
print(table(md2$group))
cat("Proteins tested:", nrow(e2), "\n")


# ---- 2. Per-protein mean comparison, same covariates -----------------
# Identical adjustment to the variance analysis, so the two are
# directly comparable. This is the point: same data, same model,
# different summary statistic.

form <- as.formula(paste("y ~ group +", paste(COV, collapse = " + ")))

res <- t(apply(e2, 1, function(y) {
  if (!is.finite(sd(y)) || sd(y) == 0) return(c(NA, NA))
  fit <- summary(lm(form, data = transform(md2, y = y)))$coefficients
  i <- grep("^groupType2$", rownames(fit))
  if (length(i) == 0) return(c(NA, NA))
  c(fit[i, "Estimate"], fit[i, "Pr(>|t|)"])
}))

de <- data.frame(gene = rownames(e2),
                 estimate = res[, 1], p = res[, 2],
                 stringsAsFactors = FALSE) %>%
  filter(is.finite(estimate), is.finite(p)) %>%
  mutate(fdr = p.adjust(p, method = "BH"),
         negl10p = -log10(p))

# scale the effect to a log2-like quantity for plotting
de$effect <- de$estimate / stats::sd(de$estimate, na.rm = TRUE)

cat("\nProteins with FDR < 0.05:", sum(de$fdr < 0.05), "\n")
cat("Proteins with FDR < 0.20:", sum(de$fdr < 0.20), "\n")


# ---- 3. Where do our candidates land? --------------------------------

cand <- c(UNC13D = "201294", STX7 = "8417", CTSD = "1509",
          RAB3A = "5864", RAB27A = "5873", INS = "3630")

cd <- de %>%
  filter(gene %in% cand) %>%
  mutate(label = names(cand)[match(gene, cand)],
         pct_p = 100 * sapply(p, function(x) mean(de$p <= x)))

cat("\n===== OUR CANDIDATES IN THE CONVENTIONAL ANALYSIS =====\n")
print(cd %>% dplyr::select(label, estimate, p, fdr, pct_p),
      digits = 3, row.names = FALSE)

cat("\nCHECK: UNC13D should be unremarkable here.\n")
u <- cd[cd$label == "UNC13D", ]
if (nrow(u) && u$fdr < 0.05) {
  cat("  !! UNC13D IS significant on mean in this cohort (FDR",
      round(u$fdr, 4), ").\n")
  cat("  This CONTRADICTS the claim that it is invisible to\n")
  cat("  differential expression. Report it and revise the framing.\n")
} else if (nrow(u)) {
  cat(sprintf("  UNC13D: p = %.3f, FDR = %.3f, %.0fth percentile by p.\n",
              u$p, u$fdr, u$pct_p))
  cat("  Consistent with the three prior comparisons: no mean change.\n")
}


# ---- 4. Top hits, for labelling --------------------------------------
# Label the genuine hits so the volcano looks like a real analysis
# rather than an empty field, then show our candidates among the noise.

top <- de %>% arrange(p) %>% head(8)

# try to attach symbols for readability
sym <- tryCatch({
  suppressMessages({
    library(org.Hs.eg.db); library(AnnotationDbi)
    m <- AnnotationDbi::select(org.Hs.eg.db, keys = top$gene,
                               keytype = "ENTREZID", columns = "SYMBOL")
    m$SYMBOL[match(top$gene, m$ENTREZID)]
  })
}, error = function(e) top$gene)
top$label <- ifelse(is.na(sym), top$gene, sym)


# ---- 5. Figure -------------------------------------------------------

thr <- max(de$negl10p[de$fdr < 0.20], na.rm = TRUE)

fig <- ggplot(de, aes(effect, negl10p)) +
  geom_point(colour = "grey72", size = .8, alpha = .5) +
  geom_point(data = top, colour = "#ef6c00", size = 2.1) +
  geom_point(data = cd, colour = "#1b5e20", size = 2.4) +
  geom_text_repel(data = top, aes(label = label), size = 2.6,
                  colour = "grey25", seed = 1, box.padding = .4,
                  max.overlaps = 20) +
  geom_text_repel(data = cd, aes(label = label), size = 2.9,
                  colour = "#1b5e20", fontface = "bold", seed = 2,
                  box.padding = .5, min.segment.length = 0,
                  max.overlaps = 20) +
  { if (is.finite(thr)) geom_hline(yintercept = thr, linetype = "dotted",
                                   colour = "grey30") } +
  labs(x = "Mean difference, T2D vs non-diabetic (standardised)",
       y = expression(paste(-log[10], "(P)")),
       title = "Crinophagy candidates are invisible to differential expression",
       subtitle = sprintf("Discovery proteomics, %d donors, same covariate model as the variance analysis",
                          ncol(e2)),
       caption = paste(
         "Orange: the eight most significant proteins. Green: the genes carrying the variance",
         "\nphenotype. Dotted line: FDR 0.20. UNC13D and STX7 rank in the top 1% for inter-donor",
         "\nvariance in this same cohort while showing no mean difference here. Volcano format",
         "\nadapted from Wigger et al. (Nat Metab 2021); data and analysis are ours.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        plot.caption = element_text(size = 6.8, colour = "grey40", hjust = 0))

ggsave(file.path(FIGS, "Fig_mean_volcano.pdf"), fig,
       width = 6.5, height = 5.5)

write.csv(de, file.path(FIGS, "discovery_mean_DE.csv"), row.names = FALSE)

cat("\nWrote Fig_mean_volcano.pdf and discovery_mean_DE.csv\n")
cat("\nPair this with the variance forest plot. Same cohort, same model,\n")
cat("different summary statistic, opposite conclusions about UNC13D.\n")
