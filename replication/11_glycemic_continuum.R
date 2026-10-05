# =====================================================================
# UNC13D variance across the full glycemic continuum
#
# WHERE THIS SITS. Established so far:
#   UNC13D  proteomics SD ratio 2.065 (discovery, n=123)
#           RNA-seq    SD ratio 1.918 (replication, ND vs T2D)
#           99.0th percentile genome-wide
#           +batch          1.918
#           +batch+immune   1.849
#           +batch+immune+beta 1.870   <- stable across all three
#   RAB3A   1.725 discovery -> 1.089 replication (does NOT replicate)
#   RAB27A  5.2nd percentile: variance REDUCED in T2D
#   GABARAP 99.8th percentile, and it is a comparator-set member
#
# THE QUESTION HERE. GSE164416 has ND, IFG, IGT, T3cD and T2D. Two
# things follow:
#
#   1. IGT tells us WHEN heterogeneity emerges. If IGT is already
#      elevated, this precedes overt diabetes, which is a better
#      claim than an end-stage association.
#
#   2. T3cD is the disease-specificity control. It is pancreatogenic
#      diabetes: diabetes secondary to pancreatic disease, NOT T2D.
#      If T3cD variance is as high as T2D, this is general islet
#      injury and the paper's framing must change. That is the
#      result that could overturn the central claim, which is why
#      this analysis runs before any writing.
#
# THE OBSTACLE. Batch is confounded with diagnosis (chi-sq p=0.00068).
# The 24/11/2016 batch has 13 IGT and 14 T3cD but ZERO ND. A
# within-batch comparison against ND will therefore drop donors from
# those groups. The script reports surviving n per group BEFORE any
# ratio, and refuses to fit a group below MIN_N.
#
# Requires: mat, pd2, cnt, pd from scripts 09/10 in the session.
# =====================================================================

library(DESeq2)
library(dplyr)
library(car)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd("/Users/mmani3/Downloads/Human Islets")
REP <- "replication_GSE164416"

UNC13D  <- "ENSG00000092929"
GABARAP <- "ENSG00000170296"
RAB3A   <- "ENSG00000105649"
RAB27A  <- "ENSG00000069974"

MIN_N   <- 8      # below this we report "underpowered", we do not fit
NBOOT   <- 5000
set.seed(1)


# ---- 1. Rebuild the full cohort (not just ND/T2D) -------------------

cnt_all <- read.delim(gzfile(file.path(REP, "GSE164416",
                            "GSE164416_DP_htseq_counts.txt.gz")),
                      row.names = 1, check.names = FALSE)

pd_all <- read.csv(file.path(REP, "gse164416_metadata.csv"),
                   stringsAsFactors = FALSE)
pd_all$dp <- sub(".*_(DP[0-9]+)_.*", "\\1", pd_all$title)
names(pd_all)[names(pd_all) == "diabetes.status"] <- "group"
pd_all <- pd_all[match(colnames(cnt_all), pd_all$dp), ]
rownames(pd_all) <- pd_all$dp
pd_all$batch <- factor(pd_all$library_prep_date)

# Drop IFG (n=3) and the 2-sample 2013 batch: neither can be adjusted
keep <- pd_all$group %in% c("ND", "IGT", "T3cD", "T2D") &
        pd_all$batch != "26/02/2013"
pd_all <- droplevels(pd_all[keep, ])
cnt_all <- cnt_all[, rownames(pd_all)]

pd_all$group <- factor(pd_all$group, levels = c("ND", "IGT", "T3cD", "T2D"))

cat("===== COHORT BEFORE BATCH FILTERING =====\n")
print(table(pd_all$group, pd_all$batch))


# ---- 2. VST on the full cohort --------------------------------------
# blind = TRUE: the transformation must not borrow information from
# the group labels, or it could partly encode the effect being tested.

cnt_all <- cnt_all[rowSums(cnt_all >= 10) >= 10, ]
dds <- DESeqDataSetFromMatrix(round(as.matrix(cnt_all)), pd_all,
                              ~ batch + group)
vsd <- vst(dds, blind = TRUE)
m_all <- assay(vsd)
cat("\nVST matrix:", dim(m_all), "\n")


# ---- 3. Covariates: immune infiltration and beta-cell content -------
# Same construction as the ND/T2D analysis, recomputed on this matrix
# so the two are directly comparable.

imm_sym  <- c("MS4A1","PTPRC","CD3E","CD79A","CD79B","MZB1","JCHAIN",
              "CCL19","CCL18","CMA1","CTSG","TPSAB1","CPA3")
beta_sym <- c("INS","IAPP","NKX6-1","PDX1","MAFA","G6PC2","ADCYAP1",
              "SLC30A8")

ens_of <- function(syms) {
  a <- AnnotationDbi::select(org.Hs.eg.db, keys = syms,
                             keytype = "SYMBOL", columns = "ENSEMBL")
  unique(a$ENSEMBL[!is.na(a$ENSEMBL) & a$ENSEMBL %in% rownames(m_all)])
}

# immunoglobulin genes dominated the high-variance tail, so include them
ig_ens <- rownames(m_all)[
  rownames(m_all) %in% ens_of(grep("^IG[HKL]",
    keys(org.Hs.eg.db, "SYMBOL"), value = TRUE))]

imm_ens  <- unique(c(ens_of(imm_sym), ig_ens))
beta_ens <- ens_of(beta_sym)

cat("Immune score genes:", length(imm_ens),
    "  Beta score genes:", length(beta_ens), "\n")

pd_all$imm  <- colMeans(m_all[imm_ens,  rownames(pd_all), drop = FALSE])
pd_all$beta <- colMeans(m_all[beta_ens, rownames(pd_all), drop = FALSE])

cat("\nCovariate SD by group (both are MORE variable in disease,\n")
cat("which is exactly why they must be adjusted for):\n")
print(rbind(immune = tapply(pd_all$imm,  pd_all$group, sd),
            beta   = tapply(pd_all$beta, pd_all$group, sd)))


# ---- 4. Pairwise comparisons against ND -----------------------------
# Each disease group is compared with ND separately, restricted to
# batches containing both. That is why n differs between comparisons:
# it is a consequence of the confounded design, not a choice.

compare_to_ND <- function(ens, label, grp) {
  sub <- pd_all[pd_all$group %in% c("ND", grp), ]
  sub <- droplevels(sub)

  tab <- table(sub$batch, sub$group)
  ok  <- rownames(tab)[rowSums(tab > 0) == 2]
  sub <- droplevels(sub[sub$batch %in% ok, ])

  n_nd <- sum(sub$group == "ND")
  n_gr <- sum(sub$group == grp)

  base <- data.frame(gene = label, comparison = paste(grp, "vs ND"),
                     n_ND = n_nd, n_grp = n_gr,
                     sd_ratio = NA, ci_lo = NA, ci_hi = NA,
                     p = NA, note = "", stringsAsFactors = FALSE)

  if (n_nd < MIN_N || n_gr < MIN_N) {
    base$note <- sprintf("UNDERPOWERED (min n = %d)", MIN_N)
    return(base)
  }

  sub$y <- m_all[ens, rownames(sub)]
  # full model: batch + immune + beta-cell content
  sub$r <- residuals(lm(y ~ batch + imm + beta, data = sub))

  a <- sub$r[sub$group == grp]
  b <- sub$r[sub$group == "ND"]

  boot <- replicate(NBOOT, {
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE))
  })

  base$sd_ratio <- sd(a) / sd(b)
  base$ci_lo <- quantile(boot, .025)
  base$ci_hi <- quantile(boot, .975)
  base$p <- leveneTest(r ~ group, data = sub)$`Pr(>F)`[1]
  base$note <- if (base$ci_lo > 1) "CI excludes 1" else ""
  base
}

genes <- c(UNC13D = UNC13D, GABARAP = GABARAP,
           RAB3A = RAB3A, RAB27A = RAB27A)

res <- do.call(rbind, lapply(names(genes), function(g)
  do.call(rbind, lapply(c("IGT", "T3cD", "T2D"),
                        function(x) compare_to_ND(genes[[g]], g, x)))))

cat("\n===== VARIANCE RATIO vs ND (batch + immune + beta adjusted) =====\n")
print(res, digits = 4, row.names = FALSE)

write.csv(res, file.path(REP, "glycemic_continuum_variance.csv"),
          row.names = FALSE)


# ---- 5. Read the pattern --------------------------------------------

u <- res[res$gene == "UNC13D", ]
cat("\n===== INTERPRETATION =====\n")

get <- function(cmp) {
  v <- u$sd_ratio[u$comparison == cmp]
  if (length(v) == 0 || is.na(v)) NA else v
}
r_igt  <- get("IGT vs ND")
r_t3cd <- get("T3cD vs ND")
r_t2d  <- get("T2D vs ND")

cat(sprintf("  IGT  %s\n", ifelse(is.na(r_igt),  "underpowered",
                                  sprintf("%.3f", r_igt))))
cat(sprintf("  T3cD %s\n", ifelse(is.na(r_t3cd), "underpowered",
                                  sprintf("%.3f", r_t3cd))))
cat(sprintf("  T2D  %s\n", ifelse(is.na(r_t2d),  "underpowered",
                                  sprintf("%.3f", r_t2d))))

if (!is.na(r_t3cd) && !is.na(r_t2d)) {
  if (r_t3cd >= r_t2d * 0.85) {
    cat("\n  T3cD is as elevated as T2D.\n")
    cat("  => This looks like GENERAL ISLET INJURY, not T2D-specific\n")
    cat("     biology. The paper's framing must change accordingly.\n")
  } else if (r_t3cd < 1.3) {
    cat("\n  T3cD is near baseline while T2D is elevated.\n")
    cat("  => Supports T2D SPECIFICITY. This is the strong outcome.\n")
  } else {
    cat("\n  T3cD is intermediate. Report it descriptively; do not\n")
    cat("     force a specificity claim the data do not carry.\n")
  }
}

if (!is.na(r_igt)) {
  if (r_igt > 1.3) {
    cat("\n  IGT already elevated => heterogeneity emerges BEFORE overt\n")
    cat("     diabetes. Materially more interesting than an end-stage\n")
    cat("     association.\n")
  } else {
    cat("\n  IGT near baseline => a late feature, appearing with overt\n")
    cat("     T2D rather than during glucose intolerance.\n")
  }
}

cat("\n  GABARAP is shown alongside deliberately. If it tracks UNC13D\n")
cat("  across every group, the comparator problem persists and the\n")
cat("  crinophagy-specific claim cannot be made.\n")


# ---- 6. Figure ------------------------------------------------------

pl <- res %>%
  filter(!is.na(sd_ratio)) %>%
  mutate(comparison = factor(comparison,
                             levels = c("IGT vs ND", "T3cD vs ND",
                                        "T2D vs ND")))

p <- ggplot(pl, aes(x = comparison, y = sd_ratio, colour = gene)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(ymin = ci_lo, ymax = ci_hi),
                  position = position_dodge(width = .5)) +
  labs(x = NULL, y = "SD ratio vs ND (95% bootstrap CI)",
       title = "Variance across the glycemic continuum",
       subtitle = "Adjusted for batch, immune infiltration, beta-cell content",
       colour = NULL) +
  theme_minimal(base_size = 12)

ggsave(file.path(REP, "fig_glycemic_continuum.pdf"), p,
       width = 7, height = 5)

cat("\nWrote glycemic_continuum_variance.csv and fig_glycemic_continuum.pdf\n")
