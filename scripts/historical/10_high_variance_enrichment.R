# =====================================================================
# What company does UNC13D keep?
#
# Genome-wide T2D/ND variance ratios (GSE164416, batch-adjusted VST):
#   GABARAP  2.243  99.8th percentile   <- a "negative control"
#   UNC13D   1.918  99.0th percentile   <- replicates the proteomics
#   RAB3A    1.089  49.8th             <- does NOT replicate
#   RAB27A   0.762   5.2nd             <- LOWER variance in T2D
#
# UNC13D being top 1% is meaningful, but ~150 genes sit above it and
# one of them is a comparator-set member. This script asks what those
# genes are. Three outcomes, each pointing somewhere different:
#
#   (a) They cluster into a coherent biology (stress, immune, ER,
#       secretory). Then UNC13D sits inside a real process and GABARAP
#       may belong there too, which explains the comparator failure
#       without special pleading.
#   (b) They are dominated by technical signatures (mitochondrial,
#       ribosomal, hemoglobin, low-expression genes). Then residual
#       batch or contamination is inflating variance and neither
#       result survives.
#   (c) No structure at all. Then the tail is noise and UNC13D's
#       percentile is less impressive than it looks.
#
# Requires: mat, pd2, ratio_all in the session (from script 09 + the
# genome-wide block). Re-derives ratio_all if absent.
# =====================================================================

library(dplyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))
REP <- "replication_GSE164416"

# ---- 0. Recompute ratios if needed ----------------------------------

if (!exists("ratio_all")) {
  ratio_all <- apply(mat, 1, function(y) {
    r <- residuals(lm(y ~ pd2$batch))
    sd(r[pd2$group == "T2D"]) / sd(r[pd2$group == "ND"])
  })
  ratio_all <- ratio_all[is.finite(ratio_all)]
}

cat("Genes with a finite ratio:", length(ratio_all), "\n")
cat("Median ratio:", round(median(ratio_all), 3), "\n\n")


# ---- 1. Define the high-variance set --------------------------------
# Everything at or above UNC13D's ratio.

unc <- "ENSG00000092929"
thr <- ratio_all[unc]
top <- names(ratio_all)[ratio_all >= thr]

cat("Threshold (UNC13D ratio):", round(thr, 3), "\n")
cat("Genes at or above it:", length(top), "\n")

ann <- AnnotationDbi::select(org.Hs.eg.db, keys = top,
                             keytype = "ENSEMBL",
                             columns = c("SYMBOL", "ENTREZID", "GENENAME"))
ann <- ann[!duplicated(ann$ENSEMBL), ]
ann$ratio <- ratio_all[ann$ENSEMBL]
ann <- ann[order(-ann$ratio), ]

cat("\nTop 40 by variance ratio:\n")
print(head(ann[, c("SYMBOL", "ratio")], 40), row.names = FALSE)

write.csv(ann, file.path(REP, "high_variance_genes.csv"), row.names = FALSE)


# ---- 2. TECHNICAL ARTEFACT CHECK (do this before enrichment) --------
# If the tail is mitochondrial, ribosomal or low-expression genes,
# it is a technical signature and no amount of enrichment analysis
# will make it biology.

sym <- ann$SYMBOL[!is.na(ann$SYMBOL)]

classes <- list(
  mitochondrial = grepl("^MT-", sym),
  ribosomal     = grepl("^RP[SL]", sym),
  hemoglobin    = grepl("^HB[ABDGQZ]", sym),
  immunoglobulin= grepl("^IG[HKL]", sym),
  pancreatic_exo= sym %in% c("PRSS1","PRSS2","CPA1","CPA2","CPB1","CELA2A",
                             "CELA3A","CTRB1","CTRB2","CLPS","PNLIP","AMY2A",
                             "AMY2B","SYCN","REG1A","REG1B","REG3A")
)

cat("\n===== TECHNICAL CLASSES IN THE HIGH-VARIANCE SET =====\n")
for (nm in names(classes)) {
  n <- sum(classes[[nm]])
  cat(sprintf("  %-16s %3d of %d (%.1f%%)\n",
              nm, n, length(sym), 100 * n / length(sym)))
}
cat("\n  Exocrine contamination matters most here: LCM islets from\n")
cat("  surgical specimens vary in acinar carryover, which inflates\n")
cat("  variance for exocrine genes and can spill over broadly.\n")

# Is high variance simply tracking low expression?
mean_expr <- rowMeans(mat)
common <- intersect(names(ratio_all), names(mean_expr))
ct <- cor.test(mean_expr[common], ratio_all[common], method = "spearman")
cat(sprintf("\n  Spearman(mean expression, variance ratio) = %.3f, p = %.2g\n",
            ct$estimate, ct$p.value))
cat("  A strong negative correlation would mean the tail is just\n")
cat("  low-expressed genes, and the percentile is uninformative.\n")

# Where do our genes sit on expression?
for (g in c(UNC13D = unc, GABARAP = "ENSG00000170296",
            RAB3A = "ENSG00000105649", RAB27A = "ENSG00000069974")) {
  if (g %in% names(mean_expr)) {
    pc <- 100 * mean(mean_expr <= mean_expr[g])
    cat(sprintf("  %-8s mean expr percentile = %.1f\n",
                names(which(c(UNC13D = unc, GABARAP = "ENSG00000170296",
                              RAB3A = "ENSG00000105649",
                              RAB27A = "ENSG00000069974") == g))[1], pc))
  }
}


# ---- 3. Enrichment --------------------------------------------------
# Background = all tested genes, NOT the whole genome. Using the
# genome as background here would produce false enrichment for
# anything islet-expressed.

run_enrich <- function() {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    cat("\n  clusterProfiler not installed. Install with:\n")
    cat('    BiocManager::install(c("clusterProfiler","org.Hs.eg.db"))\n')
    return(NULL)
  }
  library(clusterProfiler)

  bg <- AnnotationDbi::select(org.Hs.eg.db, keys = names(ratio_all),
                              keytype = "ENSEMBL", columns = "ENTREZID")
  bg <- unique(bg$ENTREZID[!is.na(bg$ENTREZID)])
  fg <- unique(ann$ENTREZID[!is.na(ann$ENTREZID)])

  cat("\n===== ENRICHMENT =====\n")
  cat("Foreground:", length(fg), "  Background:", length(bg), "\n")

  go <- enrichGO(fg, OrgDb = org.Hs.eg.db, universe = bg,
                 ont = "BP", pAdjustMethod = "BH",
                 pvalueCutoff = 0.05, qvalueCutoff = 0.2,
                 readable = TRUE)

  if (is.null(go) || nrow(as.data.frame(go)) == 0) {
    cat("\n  No GO BP terms enriched. That is outcome (c): the tail\n")
    cat("  has no coherent structure, and UNC13D's percentile should\n")
    cat("  be reported without a mechanistic story attached.\n")
  } else {
    df <- as.data.frame(go)
    cat("\nTop 20 GO BP terms:\n")
    print(head(df[, c("Description", "GeneRatio", "p.adjust", "Count")], 20),
          row.names = FALSE)
    write.csv(df, file.path(REP, "high_variance_GO.csv"), row.names = FALSE)

    p <- ggplot(head(df, 15),
                aes(x = Count, y = reorder(Description, Count),
                    fill = -log10(p.adjust))) +
      geom_col() +
      scale_fill_gradient(low = "#a5d6a7", high = "#1b5e20") +
      labs(x = "Genes", y = NULL,
           title = "Processes over-represented among high-variance genes",
           fill = "-log10 q") +
      theme_minimal(base_size = 11)
    ggsave(file.path(REP, "fig_high_variance_GO.pdf"), p,
           width = 8, height = 6)
  }
  invisible(go)
}

run_enrich()


# ---- 4. Is GABARAP's set behaving as a set? -------------------------
# If the whole ATG8 family sits high, the comparator set was simply
# a poor choice for a variance contrast, which is an honest thing to
# report. If only GABARAP is high, it is a single outlier.

atg8 <- c(GABARAP = "ENSG00000170296", GABARAPL1 = "ENSG00000139112",
          GABARAPL2 = "ENSG00000034713", MAP1LC3A = "ENSG00000101460",
          MAP1LC3B = "ENSG00000140941")

cat("\n===== ATG8 FAMILY PERCENTILES =====\n")
for (nm in names(atg8)) {
  e <- atg8[[nm]]
  if (e %in% names(ratio_all)) {
    cat(sprintf("  %-10s ratio=%.3f  percentile=%.1f\n",
                nm, ratio_all[e], 100 * mean(ratio_all <= ratio_all[e])))
  } else {
    cat(sprintf("  %-10s not in matrix\n", nm))
  }
}


# ---- 5. Honest summary ----------------------------------------------

cat("\n===== WHAT TO CONCLUDE =====\n")
cat("
  Established so far:
    UNC13D  variance elevated in T2D, replicated across proteomics
            (2.065) and transcriptomics (1.918), 99th percentile
            genome-wide. The MS dropout confound does not apply to
            RNA-seq, so the replication answers it.
    RAB3A   discovery result (1.725, p=0.005) does NOT replicate
            (1.089, 50th percentile). Report both.
    RAB27A  5th percentile: variance is REDUCED in T2D. Exploratory,
            not pre-specified, but argues against generalised
            destabilisation of secretory machinery.

  Unresolved:
    GABARAP exceeds UNC13D and is a comparator-set member, so
    crinophagy-specific variance is NOT established. Sections 2-4
    above decide whether that is technical, whether the ATG8 family
    was a poor comparator choice, or whether the tail is noise.
    Revising the comparator after seeing results is defensible only
    if stated openly as a post-hoc change.
")
