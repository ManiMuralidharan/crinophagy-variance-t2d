# =====================================================================
# Replication: RAB3A variance heterogeneity in GSE164416
#
# Discovery (HumanIslets proteomics, n=123 after adjustment):
#   RAB3A    SD ratio T2D/ND = 1.725, 95% CI [1.088, 2.261]
#            Brown-Forsythe p = 0.0051, permutation p = 0.0047
#            0/10 macroautophagy controls reached p<0.05
#   UNC13D   SD ratio 2.065, but detection confounded (chi-sq p=0.015)
#
# This cohort (Wigger/Solimena, Nat Metab 2021): LCM islets from
# pancreatectomised living donors. Different tissue handling, different
# platform, different donors. The paper's own headline is that diabetic
# islets are heterogeneous where non-diabetic ones are not, so the
# phenomenon has independent precedent.
#
# THREE THINGS THIS COHORT FORCES:
#   1. Batch is confounded with diagnosis (chi-sq p = 0.00068).
#      library_prep_date MUST be adjusted for. The 2016 batch has
#      13 IGT and 14 T3cD but zero ND.
#   2. Raw HTSeq counts have a built-in mean-variance relationship.
#      VST is required; CPM or log2 counts would manufacture the
#      very effect we are testing.
#   3. T3cD is pancreatogenic diabetes, mechanistically distinct from
#      T2D. Excluded from the primary comparison, kept as exploratory.
#      IFG (n=3) dropped.
#
# NOTE: transcript variance and protein variance are different
# quantities. A null here does not refute the protein finding.
# =====================================================================

library(DESeq2)
library(dplyr)
library(car)
library(ggplot2)

setwd("/Users/mmani3/Downloads/Human Islets")
REP <- "replication_GSE164416"

# ---- 1. Load ---------------------------------------------------------

cnt <- read.delim(gzfile(file.path(REP, "GSE164416",
                                   "GSE164416_DP_htseq_counts.txt.gz")),
                  row.names = 1, check.names = FALSE)

pd <- read.csv(file.path(REP, "gse164416_metadata.csv"),
               stringsAsFactors = FALSE)

# Extract DP id from titles like "Islets_DP002_T3cD"
pd$dp <- sub(".*_(DP[0-9]+)_.*", "\\1", pd$title)
names(pd)[names(pd) == "diabetes.status"] <- "group"
if (!"group" %in% names(pd)) names(pd)[1] <- "group"

stopifnot(all(colnames(cnt) %in% pd$dp))
pd <- pd[match(colnames(cnt), pd$dp), ]
rownames(pd) <- pd$dp

pd$group <- factor(pd$group)
pd$batch <- factor(pd$library_prep_date)

cat("Samples:", ncol(cnt), "  Genes:", nrow(cnt), "\n")
print(table(pd$group, pd$batch))


# ---- 2. Restrict the cohort -----------------------------------------
# ND vs T2D only. Drop the 2013 batch (n=2, both T2D): a two-sample
# batch cannot be adjusted for.

keep <- pd$group %in% c("ND", "T2D") &
        pd$batch != "26/02/2013"

pd2  <- droplevels(pd[keep, ])
cnt2 <- cnt[, rownames(pd2)]

cat("\nAfter restriction:\n")
print(table(pd2$group, pd2$batch))
cat("ND =", sum(pd2$group == "ND"),
    "  T2D =", sum(pd2$group == "T2D"), "\n")

# Drop batches that no longer contain both groups: they contribute
# nothing to a within-batch adjusted comparison and destabilise the fit
tab <- table(pd2$batch, pd2$group)
ok <- rownames(tab)[rowSums(tab > 0) == 2]
cat("Batches with both groups:", paste(ok, collapse = ", "), "\n")
pd2  <- droplevels(pd2[pd2$batch %in% ok, ])
cnt2 <- cnt2[, rownames(pd2)]
cat("Final: ND =", sum(pd2$group == "ND"),
    "  T2D =", sum(pd2$group == "T2D"), "\n")


# ---- 3. Filter and VST ----------------------------------------------
# VST, not CPM: it removes the mean-variance dependence that would
# otherwise masquerade as the effect of interest.

cnt2 <- cnt2[rowSums(cnt2 >= 10) >= 10, ]
cat("\nGenes after filtering:", nrow(cnt2), "\n")

dds <- DESeqDataSetFromMatrix(round(as.matrix(cnt2)), pd2, ~ batch + group)
vsd <- vst(dds, blind = TRUE)     # blind: design must not inform the
                                  # transformation used for a variance test
mat <- assay(vsd)

cat("VST matrix:", dim(mat), "\n")


# ---- 4. Gene set ----------------------------------------------------
# Ensembl IDs for the frozen sets. RAB3A is the primary; the
# macroautophagy genes are the negative-control panel that made the
# discovery result convincing.

genes <- c(
  RAB3A     = "ENSG00000105649",   # PRIMARY
  UNC13D    = "ENSG00000092929",   # secondary (detection-confounded in MS)
  UNC13B    = "ENSG00000198722",   # paralogue control (null in discovery)
  CADPS     = "ENSG00000163618",
  RAB3GAP1  = "ENSG00000115839",
  # macroautophagy negative controls
  ATG3      = "ENSG00000144848",
  ATG5      = "ENSG00000057663",
  ATG7      = "ENSG00000197548",
  BECN1     = "ENSG00000126581",
  ULK1      = "ENSG00000177169",
  WIPI2     = "ENSG00000157954",
  MAP1LC3A  = "ENSG00000101460",
  GABARAP   = "ENSG00000170296",
  GABARAPL1 = "ENSG00000139112",
  GABARAPL2 = "ENSG00000034713",
  SQSTM1    = "ENSG00000161011",
  NBR1      = "ENSG00000188554",
  CALCOCO2  = "ENSG00000136436"
)

present <- genes[genes %in% rownames(mat)]
cat("\nGenes found:", length(present), "of", length(genes), "\n")
missing <- setdiff(names(genes), names(present))
if (length(missing)) cat("Missing:", paste(missing, collapse = ", "), "\n")


# ---- 5. Variance test, batch-adjusted -------------------------------

test_gene <- function(ens, label) {
  df <- pd2
  df$y <- mat[ens, rownames(pd2)]
  df$r <- residuals(lm(y ~ batch, data = df))   # batch adjustment

  a <- df$r[df$group == "T2D"]
  b <- df$r[df$group == "ND"]

  lev <- leveneTest(r ~ group, data = df)

  set.seed(1)
  boot <- replicate(5000, {
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE))
  })

  data.frame(gene = label,
             n_T2D = length(a), n_ND = length(b),
             sd_ratio = sd(a) / sd(b),
             ci_lo = quantile(boot, .025),
             ci_hi = quantile(boot, .975),
             p = lev$`Pr(>F)`[1],
             row.names = NULL)
}

res <- do.call(rbind, Map(test_gene, present, names(present)))
res <- res[order(-res$sd_ratio), ]

cat("\n===== VARIANCE RATIO (T2D / ND), batch-adjusted =====\n")
print(res, digits = 4, row.names = FALSE)

cat("\n===== PRIMARY: RAB3A =====\n")
print(res[res$gene == "RAB3A", ], digits = 4, row.names = FALSE)
cat("Discovery cohort: SD ratio 1.725, 95% CI [1.088, 2.261], p = 0.0051\n")

ctrl <- res[res$gene %in% c("ATG3","ATG5","ATG7","BECN1","ULK1","WIPI2",
                            "MAP1LC3A","GABARAP","GABARAPL1","GABARAPL2",
                            "SQSTM1","NBR1","CALCOCO2"), ]
cat("\nControls reaching p<0.05:", sum(ctrl$p < 0.05), "of", nrow(ctrl),
    " (expect ~", round(0.05*nrow(ctrl), 1), "by chance)\n")
cat("Max control SD ratio:", round(max(ctrl$sd_ratio), 3), "\n")


# ---- 6. Permutation, RAB3A only -------------------------------------

df <- pd2
df$y <- mat[present["RAB3A"], rownames(pd2)]
df$r <- residuals(lm(y ~ batch, data = df))

obs <- leveneTest(r ~ group, data = df)$`F value`[1]
set.seed(1)
null <- replicate(10000, {
  d <- df; d$group <- sample(d$group)
  leveneTest(r ~ group, data = d)$`F value`[1]
})
cat("\nPermutation p (RAB3A):", (sum(null >= obs) + 1) / 10001, "\n")


# ---- 7. Sensitivity: authors' own quality subset --------------------

if ("in_ins_filtered_data_subset" %in% names(pd2)) {
  sel <- pd2$in_ins_filtered_data_subset %in% c("TRUE", TRUE)
  cat("\n===== Restricted to in_ins_filtered subset =====\n")
  cat("n =", sum(sel), "\n")
  if (sum(sel) > 20) {
    d3 <- pd2[sel, ]
    d3$y <- mat[present["RAB3A"], rownames(d3)]
    d3$r <- residuals(lm(y ~ batch, data = droplevels(d3)))
    a <- d3$r[d3$group == "T2D"]; b <- d3$r[d3$group == "ND"]
    cat("ND =", length(b), " T2D =", length(a),
        " SD ratio =", round(sd(a)/sd(b), 3),
        " p =", round(leveneTest(r ~ group, data = d3)$`Pr(>F)`[1], 4), "\n")
  }
}


# ---- 8. Exploratory: T3cD -------------------------------------------
# If RAB3A variance rises in T2D but not T3cD, that is a specificity
# result worth having. Interpret cautiously: different disease, and
# batch structure differs again.

keep3 <- pd$group %in% c("ND", "T3cD") & pd$batch != "26/02/2013"
pd3 <- droplevels(pd[keep3, ])
tab3 <- table(pd3$batch, pd3$group)
ok3 <- rownames(tab3)[rowSums(tab3 > 0) == 2]
pd3 <- droplevels(pd3[pd3$batch %in% ok3, ])

if (nrow(pd3) > 20) {
  c3 <- cnt[rownames(mat), rownames(pd3)]
  d3 <- DESeqDataSetFromMatrix(round(as.matrix(c3)), pd3, ~ batch + group)
  m3 <- assay(vst(d3, blind = TRUE))
  pd3$y <- m3[present["RAB3A"], rownames(pd3)]
  pd3$r <- residuals(lm(y ~ batch, data = pd3))
  a <- pd3$r[pd3$group == "T3cD"]; b <- pd3$r[pd3$group == "ND"]
  cat("\n===== EXPLORATORY: T3cD vs ND =====\n")
  cat("ND =", length(b), " T3cD =", length(a),
      " SD ratio =", round(sd(a)/sd(b), 3),
      " p =", round(leveneTest(r ~ group, data = pd3)$`Pr(>F)`[1], 4), "\n")
}


# ---- 9. Save --------------------------------------------------------

write.csv(res, file.path(REP, "replication_variance_results.csv"),
          row.names = FALSE)

df$group <- factor(df$group, levels = c("ND", "T2D"))
p <- ggplot(df, aes(group, r, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .6) +
  geom_jitter(width = .15, alpha = .6) +
  labs(x = NULL, y = "RAB3A (VST, batch-adjusted residual)",
       title = "RAB3A dispersion, replication cohort") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(REP, "fig_rab3a_replication.pdf"), p, width = 5, height = 5)

cat("\nWrote replication_variance_results.csv and fig_rab3a_replication.pdf\n")
