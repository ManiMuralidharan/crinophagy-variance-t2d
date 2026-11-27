# =====================================================================
# 32. Sensitivity analysis: beta-cell content score without INS
#
# WHY. INS is a member of the pre-registered 31-gene crinophagy set and
# also one of the 8 genes in the beta-cell content score used as a
# covariate in Model 2. INS is therefore partly regressed on itself,
# and its attenuation after adjustment (1.159 -> 0.974) cannot be
# distinguished from that circularity. Methods 4.6 currently discloses
# this as a limitation. This script removes it.
#
# WHAT IT DOES. Rebuilds the beta-cell score from 7 genes instead of 8,
# refits Model 2, recomputes the genome-wide reference distribution
# under the new model, and reports every key gene old versus new.
#
# WHAT WOULD CHANGE THE PAPER. If UNC13D and STX7 hold their ratios and
# percentiles, nothing changes except that the limitation paragraph is
# replaced by a result. If INS moves back toward 1.159, you report that
# its adjusted value was partly self-adjusted, which is what you already
# suspect. Either outcome is reportable.
#
# RUN:  source("32_ins_sensitivity.R")
# =====================================================================

library(car)

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))

# ---- 0. Session objects ---------------------------------------------
if (!all(sapply(c("mat", "pd2", "bg2"), exists))) {
  load("session_gse164416.RData")
}
stopifnot(exists("mat"), exists("pd2"), exists("bg2"))
stopifnot(all(c("batch", "imm", "beta", "group") %in% names(pd2)))

set.seed(1)
samples <- rownames(pd2)
cat("Donors:", nrow(pd2), " Genes:", nrow(mat), "\n")
print(table(pd2$group))


# ---- 1. Resolve the 7 remaining beta-cell marker genes ---------------
# Resolved programmatically. Do not hand-type Ensembl accessions.

BETA8 <- c("INS", "IAPP", "NKX6-1", "PDX1", "MAFA", "G6PC2",
           "ADCYAP1", "SLC30A8")
BETA7 <- setdiff(BETA8, "INS")

suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi)
})
mp <- AnnotationDbi::select(org.Hs.eg.db, keys = BETA7,
                            keytype = "SYMBOL", columns = "ENSEMBL")
detach("package:org.Hs.eg.db", unload = TRUE)
detach("package:AnnotationDbi", unload = TRUE)

mp  <- mp[!is.na(mp$ENSEMBL) & mp$ENSEMBL %in% rownames(mat), ]
ids <- unique(mp$ENSEMBL)

cat("\nBeta-cell markers resolved and present:",
    length(unique(mp$SYMBOL)), "of", length(BETA7), "\n")
print(mp)
missing <- setdiff(BETA7, mp$SYMBOL)
if (length(missing))
  cat("NOT FOUND (report this):", paste(missing, collapse = ", "), "\n")
stopifnot(length(ids) >= 5)


# ---- 2. Rebuild the score -------------------------------------------
# Same construction as the original: unweighted mean of marker
# expression on the variance-stabilised scale.

pd2$beta7 <- colMeans(mat[ids, samples, drop = FALSE])

cat("\nOld (8-gene) vs new (7-gene) beta score:\n")
cat("  Pearson  r =", round(cor(pd2$beta, pd2$beta7), 4), "\n")
cat("  Spearman r =",
    round(cor(pd2$beta, pd2$beta7, method = "spearman"), 4), "\n")
cat("  SD by group (7-gene):\n")
print(tapply(pd2$beta7, pd2$group, sd))


# ---- 3. Genome-wide reference distribution under the new model ------
# Recomputed, because a percentile is only meaningful against a
# background built from the same model.

bg_under <- function(formula) {
  X   <- model.matrix(formula, data = pd2)
  Y   <- t(mat[, samples, drop = FALSE])          # samples x genes
  res <- lm.fit(X, Y)$residuals
  a <- res[pd2$group == "T2D", , drop = FALSE]
  b <- res[pd2$group == "ND",  , drop = FALSE]
  apply(a, 2, sd) / apply(b, 2, sd)
}

bg2b <- bg_under(~ batch + imm + beta7)

cat("\n===== GENOME-WIDE BACKGROUND =====\n")
cat(sprintf("  Model 2  (8-gene score): median %.4f, 99th pct %.4f\n",
            median(bg2), quantile(bg2, .99)))
cat(sprintf("  Model 2b (7-gene score): median %.4f, 99th pct %.4f\n",
            median(bg2b), quantile(bg2b, .99)))
cat("  Spearman correlation of gene ranks:",
    round(cor(bg2, bg2b, method = "spearman"), 4), "\n")


# ---- 4. Per-gene comparison -----------------------------------------

one_gene <- function(ens, bg, fml, nboot = 5000) {
  d <- pd2
  d$y <- mat[ens, samples]
  d$r <- residuals(lm(update(fml, y ~ .), data = d))
  a <- d$r[d$group == "T2D"]; b <- d$r[d$group == "ND"]
  r <- sd(a) / sd(b)
  bt <- replicate(nboot,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))
  c(sd_ratio = r,
    ci_lo = unname(quantile(bt, .025)),
    ci_hi = unname(quantile(bt, .975)),
    p   = leveneTest(r ~ group, data = d)$`Pr(>F)`[1],
    pct = 100 * mean(bg <= r))
}

TARGETS <- c(UNC13D = "ENSG00000092929", STX7 = "ENSG00000079950",
             CTSD = "ENSG00000117984", GABARAP = "ENSG00000170296",
             INS = "ENSG00000254647", RAB3A = "ENSG00000105649",
             RAB27A = "ENSG00000069974")

# verify every accession against the matrix rather than trusting it
TARGETS <- TARGETS[TARGETS %in% rownames(mat)]
cat("\nTargets found:", paste(names(TARGETS), collapse = ", "), "\n")
if (length(TARGETS) < 7)
  cat("MISSING — resolve these through org.Hs.eg.db before reporting.\n")

old <- t(sapply(TARGETS, one_gene, bg = bg2,  fml = ~ batch + imm + beta))
new <- t(sapply(TARGETS, one_gene, bg = bg2b, fml = ~ batch + imm + beta7))

cmp <- data.frame(
  gene        = names(TARGETS),
  ratio_8gene = round(old[, "sd_ratio"], 3),
  ratio_7gene = round(new[, "sd_ratio"], 3),
  delta       = round(new[, "sd_ratio"] - old[, "sd_ratio"], 3),
  pct_8gene   = round(old[, "pct"], 2),
  pct_7gene   = round(new[, "pct"], 2),
  p_8gene     = signif(old[, "p"], 3),
  p_7gene     = signif(new[, "p"], 3),
  row.names   = NULL
)

cat("\n===== PER-GENE: 8-gene vs 7-gene beta score =====\n")
print(cmp, row.names = FALSE)


# ---- 5. Set-level median-shift test ---------------------------------
# Only runs if the frozen set is available. Edit the path if needed.

SETFILE <- "frozen_gene_sets_v2_clean.txt"
if (file.exists(SETFILE)) {
  gs <- read.delim(SETFILE)
  col <- grep("ENSEMBL", names(gs), ignore.case = TRUE, value = TRUE)[1]
  if (!is.na(col)) {
    crino <- intersect(as.character(gs[[col]][gs$set == "crinophagy"]),
                       names(bg2b))
    cat("\nCrinophagy genes in matrix:", length(crino), "\n")

    set_perm <- function(bg, members, nperm = 10000) {
      obs  <- median(bg[members])
      null <- replicate(nperm, median(sample(bg, length(members))))
      c(observed = obs, null_median = median(null),
        P = (sum(null >= obs) + 1) / (nperm + 1))
    }
    cat("\n===== SET-LEVEL MEDIAN SHIFT =====\n")
    cat("  8-gene score: ");
      print(round(set_perm(bg2,  crino), 4))
    cat("  7-gene score: ");
      print(round(set_perm(bg2b, crino), 4))
  }
} else {
  cat("\nSet file not found; set-level test skipped.\n")
  cat("Edit SETFILE at the top of section 5 to point at the frozen set.\n")
}


# ---- 6. Save --------------------------------------------------------

write.csv(cmp, "ins_sensitivity_comparison.csv", row.names = FALSE)
save(bg2b, file = "bg2b_beta7.RData")
cat("\nWrote ins_sensitivity_comparison.csv and bg2b_beta7.RData\n")
cat("\nREPORTING. Put the UNC13D and STX7 rows in Methods 4.6 in place\n")
cat("of the sentence saying the analysis was not performed, and the\n")
cat("full table in Table S1 or a new supplementary table.\n")
