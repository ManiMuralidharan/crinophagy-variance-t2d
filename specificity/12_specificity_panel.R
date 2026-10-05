# =====================================================================
# FINAL ANALYSIS: specificity panel
#
# This is the last analysis. After it, stop and write.
#
# What is already established:
#   UNC13D   discovery proteomics  SD ratio 2.065
#            replication RNA-seq   2.060 (batch+immune+beta adjusted)
#            99.0th percentile genome-wide
#            T2D only: IGT 1.375 (ns), T3cD 1.449 (ns), T2D 2.060 (p=0.027)
#            => T2D-specific, not general islet injury
#   RAB3A    1.725 discovery -> 1.11 replication. Does NOT replicate.
#   RAB27A   flat everywhere (0.83-1.17)
#   GABARAP  1.686 / 1.893 / 2.189 across IGT / T3cD / T2D.
#            Tracks UNC13D at every step. Survives all adjustments.
#
# THE ONE REMAINING QUESTION. GABARAP is a comparator-set member that
# behaves like UNC13D, so crinophagy-specific variance is not
# established. Is GABARAP alone, or is the whole comparator panel
# elevated? Two outcomes, two different abstracts:
#
#   (a) Only the ATG8 lipidation family (GABARAP, GABARAPL1/2,
#       MAP1LC3A) is elevated while the core machinery (ATG3/5/7,
#       BECN1, ULK1, WIPI2) and cargo receptors (SQSTM1, NBR1,
#       CALCOCO2) sit near 1.0.
#       => The comparator was not functionally independent. GABARAP
#          has documented roles at secretory granules and lysosomes
#          beyond macroautophagy, i.e. the same compartment crinophagy
#          operates in. Coherent, reportable, and it explains the
#          comparator failure without special pleading.
#
#   (b) Most of the panel is elevated.
#       => Something broad inflates variance in T2D islets and UNC13D
#          is one gene among many. Report it that way.
#
# Either way this is descriptive. The comparator set was frozen before
# analysis; we are NOT revising it to obtain a result.
#
# Requires: scripts 09-11 run in this session (m_all, pd_all).
# Rebuilds them if absent.
# =====================================================================

library(DESeq2)
library(dplyr)
library(car)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd("/Users/mmani3/Downloads/Human Islets")
REP <- "replication_GSE164416"

NBOOT <- 5000
set.seed(1)


# ---- 1. Rebuild if needed -------------------------------------------

if (!exists("m_all") || !exists("pd_all")) {
  cat("Rebuilding cohort...\n")

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

  keep <- pd_all$group %in% c("ND", "IGT", "T3cD", "T2D") &
          pd_all$batch != "26/02/2013"
  pd_all  <- droplevels(pd_all[keep, ])
  cnt_all <- cnt_all[, rownames(pd_all)]
  pd_all$group <- factor(pd_all$group,
                         levels = c("ND", "IGT", "T3cD", "T2D"))

  cnt_all <- cnt_all[rowSums(cnt_all >= 10) >= 10, ]
  dds <- DESeqDataSetFromMatrix(round(as.matrix(cnt_all)), pd_all,
                                ~ batch + group)
  m_all <- assay(vst(dds, blind = TRUE))

  ens_of <- function(syms) {
    a <- AnnotationDbi::select(org.Hs.eg.db, keys = syms,
                               keytype = "SYMBOL", columns = "ENSEMBL")
    unique(a$ENSEMBL[!is.na(a$ENSEMBL) & a$ENSEMBL %in% rownames(m_all)])
  }
  ig <- ens_of(grep("^IG[HKL]", keys(org.Hs.eg.db, "SYMBOL"), value = TRUE))
  imm_ens <- unique(c(ens_of(c("MS4A1","PTPRC","CD3E","CD79A","CD79B",
                               "MZB1","JCHAIN","CCL19","CCL18","CMA1",
                               "CTSG","TPSAB1","CPA3")), ig))
  beta_ens <- ens_of(c("INS","IAPP","NKX6-1","PDX1","MAFA","G6PC2",
                       "ADCYAP1","SLC30A8"))
  pd_all$imm  <- colMeans(m_all[imm_ens,  rownames(pd_all), drop = FALSE])
  pd_all$beta <- colMeans(m_all[beta_ens, rownames(pd_all), drop = FALSE])
}


# ---- 2. The panel, grouped by function ------------------------------
# Grouping is from the frozen v2 set plus standard pathway roles. It is
# a presentation choice, not a redefinition: every gene below was in
# the pre-registered comparator, and none is dropped.

panel <- list(
  primary = c(UNC13D = "ENSG00000092929"),

  # ATG8 family: lipidated, and several members have documented
  # non-autophagic roles at secretory granules and lysosomes
  atg8 = c(GABARAP   = "ENSG00000170296",
           GABARAPL1 = "ENSG00000139112",
           GABARAPL2 = "ENSG00000034713",
           MAP1LC3A  = "ENSG00000101460",
           MAP1LC3B  = "ENSG00000140941"),

  # core autophagy machinery: initiation and conjugation
  core = c(ATG3  = "ENSG00000144848",
           ATG5  = "ENSG00000057663",
           ATG7  = "ENSG00000197548",
           BECN1 = "ENSG00000126581",
           ULK1  = "ENSG00000177169",
           WIPI2 = "ENSG00000157954"),

  # selective autophagy cargo receptors
  receptors = c(SQSTM1   = "ENSG00000161011",
                NBR1     = "ENSG00000188554",
                CALCOCO2 = "ENSG00000136436"),

  # crinophagy set members carried through for context
  crino = c(RAB3A    = "ENSG00000105649",
            RAB27A   = "ENSG00000069974",
            CADPS    = "ENSG00000163618",
            RAB3GAP1 = "ENSG00000115839",
            UNC13B   = "ENSG00000198722")
)

flat <- unlist(panel)
names(flat) <- sub("^[a-z0-9]+\\.", "", names(flat))
grp_of <- rep(names(panel), sapply(panel, length))
names(grp_of) <- names(flat)


# ---- 3. ND vs T2D, identical adjusted model -------------------------

sub <- droplevels(pd_all[pd_all$group %in% c("ND", "T2D"), ])
tab <- table(sub$batch, sub$group)
ok  <- rownames(tab)[rowSums(tab > 0) == 2]
sub <- droplevels(sub[sub$batch %in% ok, ])

cat("\n===== COHORT =====\n")
cat("ND =", sum(sub$group == "ND"),
    "  T2D =", sum(sub$group == "T2D"), "\n")
cat("Model: ~ batch + immune + beta-cell content\n")

test_one <- function(ens, label) {
  if (!ens %in% rownames(m_all)) {
    return(data.frame(gene = label, panel = grp_of[[label]],
                      sd_ratio = NA, ci_lo = NA, ci_hi = NA,
                      p = NA, note = "not in matrix",
                      stringsAsFactors = FALSE))
  }
  d <- sub
  d$y <- m_all[ens, rownames(d)]
  d$r <- residuals(lm(y ~ batch + imm + beta, data = d))

  a <- d$r[d$group == "T2D"]
  b <- d$r[d$group == "ND"]

  boot <- replicate(NBOOT, {
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE))
  })

  data.frame(gene = label, panel = grp_of[[label]],
             sd_ratio = sd(a) / sd(b),
             ci_lo = quantile(boot, .025),
             ci_hi = quantile(boot, .975),
             p = leveneTest(r ~ group, data = d)$`Pr(>F)`[1],
             note = "", stringsAsFactors = FALSE)
}

res <- do.call(rbind, lapply(names(flat),
                             function(g) test_one(flat[[g]], g)))

# BH across the comparator panel only (the pre-registered multiple
# comparison); UNC13D is the primary test and is not part of it
is_comp <- res$panel %in% c("atg8", "core", "receptors")
res$q <- NA
res$q[is_comp] <- p.adjust(res$p[is_comp], method = "BH")

res <- res[order(-res$sd_ratio), ]

cat("\n===== FULL PANEL, ND vs T2D (fully adjusted) =====\n")
print(res, digits = 4, row.names = FALSE)

write.csv(res, file.path(REP, "specificity_panel.csv"), row.names = FALSE)


# ---- 4. Which outcome is it? ----------------------------------------

cat("\n===== WHICH OUTCOME =====\n")

med <- res %>% filter(!is.na(sd_ratio)) %>%
  group_by(panel) %>%
  summarise(n = n(), median_ratio = median(sd_ratio),
            n_sig = sum(p < 0.05, na.rm = TRUE), .groups = "drop")
print(as.data.frame(med), digits = 3, row.names = FALSE)

atg8_med <- med$median_ratio[med$panel == "atg8"]
core_med <- med$median_ratio[med$panel == "core"]
rec_med  <- med$median_ratio[med$panel == "receptors"]

if (length(atg8_med) && length(core_med)) {
  if (atg8_med > 1.4 && core_med < 1.3) {
    cat("\n  OUTCOME (a): the ATG8 lipidation family is elevated while\n")
    cat("  core machinery and cargo receptors are not.\n")
    cat("  => The comparator was not functionally independent of the\n")
    cat("     compartment under study. GABARAP-family proteins act at\n")
    cat("     secretory granules and lysosomes, the same membranes\n")
    cat("     crinophagy uses. State this plainly; it explains the\n")
    cat("     comparator failure without revising the frozen set.\n")
  } else if (mean(res$sd_ratio[is_comp] > 1.4, na.rm = TRUE) > 0.5) {
    cat("\n  OUTCOME (b): most of the comparator panel is elevated.\n")
    cat("  => Something broad inflates variance in T2D islets. UNC13D\n")
    cat("     is one gene among several. Report it that way; do not\n")
    cat("     claim specificity.\n")
  } else {
    cat("\n  MIXED: GABARAP is largely alone in the comparator set.\n")
    cat("  => Report it as an unexplained single outlier. Do not\n")
    cat("     construct a post-hoc rationale for one gene.\n")
  }
}

cat("\n  In every outcome, the UNC13D result itself stands: replicated\n")
cat("  across platforms, T2D-specific against T3cD, robust to three\n")
cat("  confounds. What changes is only how the specificity sentence\n")
cat("  in the abstract is worded.\n")


# ---- 5. Forest plot -------------------------------------------------

pl <- res %>%
  filter(!is.na(sd_ratio)) %>%
  mutate(gene = factor(gene, levels = rev(gene)),
         panel = factor(panel,
                        levels = c("primary","atg8","core",
                                   "receptors","crino")))

p <- ggplot(pl, aes(x = sd_ratio, y = gene, colour = panel)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi)) +
  scale_colour_manual(values = c(primary = "#1b5e20", atg8 = "#ef6c00",
                                 core = "#607d8b", receptors = "#90a4ae",
                                 crino = "#5c6bc0")) +
  labs(x = "SD ratio T2D / ND (95% bootstrap CI)", y = NULL,
       colour = NULL,
       title = "Specificity panel, fully adjusted",
       subtitle = "Batch + immune infiltration + beta-cell content") +
  theme_minimal(base_size = 11)

ggsave(file.path(REP, "fig_specificity_panel.pdf"), p,
       width = 7.5, height = 6)

cat("\nWrote specificity_panel.csv and fig_specificity_panel.pdf\n")
cat("\nThis completes the analysis. Next step is the draft.\n")
