# =====================================================================
# CANONICAL ANALYSIS - freeze after this
#
# Produces two clearly labelled models over one fixed 20,255-gene
# background, so every number in the manuscript comes from one version.
#
#   Model 1  ~ batch
#   Model 2  ~ batch + immune infiltration + beta-cell content
#
# WHAT IS TESTED AND WHAT IS NOT
#
#   TESTED (legitimate): the 31-gene crinophagy set was frozen with a
#   timestamped DOI before any analysis. Asking whether that set is
#   shifted relative to the genome-wide background is a pre-registered
#   test. Likewise the 14-gene comparator set.
#
#   NOT TESTED (deliberately): the "fusion arm vs cargo arm" split.
#   UNC13D and STX7 are in the fusion group BECAUSE they ranked 1st and
#   2nd. Any permutation test on a group defined after seeing the ranks
#   returns a meaningless p-value, and the frozen set contains no such
#   subdivision. The pattern is reported descriptively, labelled as
#   post-hoc, and left for independent testing.
#
#   POST-HOC COMPARATORS: INS and STX1A were added after the fact
#   because Lawlor et al. identified them as mean-shifted in T2D beta
#   cells. They are labelled as such throughout.
#
# Requires: mat, pd2, ratio_all in session (or session_gse164416.RData)
# =====================================================================

library(dplyr)
library(car)
library(ggplot2)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd("/Users/mmani3/Downloads/Human Islets")
REP <- "replication_GSE164416"

NBOOT <- 5000
NPERM <- 10000
set.seed(1)

if (!exists("mat")) load("session_gse164416.RData")


# ---- 1. Resolve every accession programmatically --------------------
# The earlier hand-typed list had VAMP7 wrong (ENSG00000124675 does not
# exist; the correct ID is ENSG00000124333). Map from symbols instead
# so this cannot recur, and record which genes are genuinely absent.

crino <- c("UNC13D","STX7","STX8","STX17","VAMP7","VAMP8","SNAP23",
           "VPS33A","VPS16","VPS11","RAB3A","RAB3GAP1","RAB27A","CADPS",
           "SYTL4","PCLO","LAMP1","LAMP2","CTSB","CTSD","CTSL",
           "ATP6V0D1","ATP6V1A","TFEB","PIKFYVE","INS","CHGA","CHGB",
           "PCSK1","PCSK2","CPE")

comp  <- c("ATG3","ATG5","ATG7","BECN1","ULK1","WIPI2","MAP1LC3A",
           "MAP1LC3B","GABARAP","GABARAPL1","GABARAPL2","SQSTM1",
           "NBR1","CALCOCO2")

posthoc <- c("STX1A")     # INS is already in the frozen crinophagy set

resolve <- function(syms) {
  m <- AnnotationDbi::select(org.Hs.eg.db, keys = syms,
                             keytype = "SYMBOL", columns = "ENSEMBL")
  m <- m[!is.na(m$ENSEMBL), ]
  m <- m[m$ENSEMBL %in% rownames(mat), ]
  m[!duplicated(m$SYMBOL), ]
}

r_crino <- resolve(crino)   %>% mutate(set = "crinophagy")
r_comp  <- resolve(comp)    %>% mutate(set = "comparator")
r_post  <- resolve(posthoc) %>% mutate(set = "post-hoc")

key <- bind_rows(r_crino, r_comp, r_post)

cat("===== ACCESSION RESOLUTION =====\n")
for (nm in c("crinophagy","comparator","post-hoc")) {
  want <- switch(nm, crinophagy = crino, comparator = comp, posthoc)
  got  <- key$SYMBOL[key$set == nm]
  miss <- setdiff(want, got)
  cat(sprintf("  %-12s %2d of %2d resolved%s\n", nm, length(got),
              length(want),
              if (length(miss)) paste0("  MISSING: ",
                                       paste(miss, collapse=", ")) else ""))
}
cat("\n  A missing gene here means it is absent from the filtered\n")
cat("  matrix, i.e. below the expression threshold. It does NOT mean\n")
cat("  the accession is wrong - all IDs are resolved from org.Hs.eg.db.\n")


# ---- 2. Covariates for Model 2 --------------------------------------

ens_of <- function(syms) {
  a <- AnnotationDbi::select(org.Hs.eg.db, keys = syms,
                             keytype = "SYMBOL", columns = "ENSEMBL")
  unique(a$ENSEMBL[!is.na(a$ENSEMBL) & a$ENSEMBL %in% rownames(mat)])
}

ig  <- ens_of(grep("^IG[HKL]", keys(org.Hs.eg.db, "SYMBOL"), value = TRUE))
imm <- unique(c(ens_of(c("MS4A1","PTPRC","CD3E","CD79A","CD79B","MZB1",
                         "JCHAIN","CCL19","CCL18","CMA1","CTSG",
                         "TPSAB1","CPA3")), ig))
bet <- ens_of(c("INS","IAPP","NKX6-1","PDX1","MAFA","G6PC2",
                "ADCYAP1","SLC30A8"))

pd2$imm  <- colMeans(mat[imm, rownames(pd2), drop = FALSE])
pd2$beta <- colMeans(mat[bet, rownames(pd2), drop = FALSE])

cat("\nImmune score:", length(imm), "genes   Beta score:",
    length(bet), "genes\n")
cat("Covariate SD by group:\n")
print(rbind(immune = tapply(pd2$imm,  pd2$group, sd),
            beta   = tapply(pd2$beta, pd2$group, sd)), digits = 3)


# ---- 3. Genome-wide backgrounds, one per model ----------------------
# Percentiles must come from the SAME model as the gene's own ratio.
# Mixing a Model-2 ratio against a Model-1 background is the error
# that made the earlier tables inconsistent.

vr <- function(y, form) {
  r <- residuals(lm(as.formula(form), data = transform(pd2, y = y)))
  s1 <- sd(r[pd2$group == "T2D"]); s0 <- sd(r[pd2$group == "ND"])
  if (!is.finite(s0) || s0 == 0) return(NA_real_)
  s1 / s0
}

cat("\nComputing genome-wide backgrounds (20k genes x 2 models)...\n")
bg1 <- apply(mat, 1, function(y) vr(y, "y ~ batch"))
bg2 <- apply(mat, 1, function(y) vr(y, "y ~ batch + imm + beta"))
bg1 <- bg1[is.finite(bg1)]
bg2 <- bg2[is.finite(bg2)]

cat("Background sizes:", length(bg1), "and", length(bg2), "\n")
cat("Median ratio  M1:", round(median(bg1),3),
    "  M2:", round(median(bg2),3), "\n")


# ---- 4. Per-gene results under both models --------------------------

one_gene <- function(ens, label, set, form, bg) {
  d <- pd2
  d$y <- mat[ens, rownames(d)]
  d$r <- residuals(lm(as.formula(form), data = d))

  a <- d$r[d$group == "T2D"]; b <- d$r[d$group == "ND"]

  boot <- replicate(NBOOT,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))

  data.frame(gene = label, set = set, ensembl = ens,
             sd_ratio = sd(a) / sd(b),
             ci_lo = quantile(boot, .025),
             ci_hi = quantile(boot, .975),
             p = leveneTest(r ~ group, data = d)$`Pr(>F)`[1],
             pct = 100 * mean(bg <= (sd(a)/sd(b))),
             row.names = NULL)
}

run_model <- function(form, bg, tag) {
  res <- do.call(rbind, lapply(seq_len(nrow(key)), function(i)
    one_gene(key$ENSEMBL[i], key$SYMBOL[i], key$set[i], form, bg)))
  # BH within the pre-registered comparator panel only; the crinophagy
  # genes are the hypothesis, the comparators are the multiple test
  res$q <- NA
  i <- res$set == "comparator"
  res$q[i] <- p.adjust(res$p[i], method = "BH")
  res$model <- tag
  res %>% arrange(desc(sd_ratio))
}

m1 <- run_model("y ~ batch",                    bg1, "M1_batch")
m2 <- run_model("y ~ batch + imm + beta",       bg2, "M2_full")

cat("\n===== MODEL 1: batch only =====\n")
print(m1 %>% select(gene, set, sd_ratio, ci_lo, ci_hi, p, pct, q),
      digits = 3, row.names = FALSE)

cat("\n===== MODEL 2: batch + immune + beta =====\n")
print(m2 %>% select(gene, set, sd_ratio, ci_lo, ci_hi, p, pct, q),
      digits = 3, row.names = FALSE)


# ---- 5. Side-by-side: what changes between models? ------------------

cmp <- m1 %>%
  select(gene, set, r1 = sd_ratio, p1 = p, pct1 = pct, q1 = q) %>%
  inner_join(m2 %>% select(gene, r2 = sd_ratio, p2 = p,
                           pct2 = pct, q2 = q), by = "gene") %>%
  mutate(delta = r2 - r1) %>%
  arrange(desc(r2))

cat("\n===== MODEL 1 vs MODEL 2 =====\n")
print(cmp, digits = 3, row.names = FALSE)

cat("\n  Genes whose ratio moves most on adjustment are those whose\n")
cat("  variance was driven by cell-composition differences rather\n")
cat("  than by cell-intrinsic variation.\n")

write.csv(cmp, file.path(REP, "canonical_both_models.csv"),
          row.names = FALSE)
write.csv(bind_rows(m1, m2), file.path(REP, "canonical_full.csv"),
          row.names = FALSE)


# ---- 6. PRE-REGISTERED SET TEST -------------------------------------
# Legitimate because both sets were frozen with a DOI before analysis.
# Question: is the set's median variance ratio shifted relative to
# random gene sets of the same size drawn from the same background?

set_test <- function(res, bg, tag) {
  cat("\n===== SET-LEVEL TEST (", tag, ") =====\n", sep = "")
  for (s in c("crinophagy", "comparator")) {
    obs_g <- res$sd_ratio[res$set == s]
    n <- length(obs_g)
    obs <- median(obs_g)
    null <- replicate(NPERM, median(sample(bg, n)))
    p_hi <- (sum(null >= obs) + 1) / (NPERM + 1)
    cat(sprintf("  %-12s n=%2d  median=%.3f  null median=%.3f  p(higher)=%.4f\n",
                s, n, obs, median(null), p_hi))
  }
  # tail enrichment: how many set members exceed the 99th percentile?
  thr <- quantile(bg, .99)
  for (s in c("crinophagy", "comparator")) {
    g <- res$sd_ratio[res$set == s]
    k <- sum(g >= thr, na.rm = TRUE)
    n <- length(g)
    pb <- binom.test(k, n, 0.01, alternative = "greater")$p.value
    cat(sprintf("  %-12s %d of %d above 99th pct (expect %.2f)  binom p=%.4g\n",
                s, k, n, 0.01 * n, pb))
  }
}

set_test(m1, bg1, "Model 1")
set_test(m2, bg2, "Model 2")


# ---- 7. DESCRIPTIVE ONLY: fusion vs cargo ---------------------------
# NO p-value is computed here, by design. UNC13D and STX7 are in the
# fusion group because they ranked 1st and 2nd; a test on a group
# defined that way is circular. Reported as an observation requiring
# independent testing.

fusion <- c("UNC13D","STX7","STX8","STX17","VAMP7","VAMP8","SNAP23",
            "VPS33A","VPS16","VPS11","CTSB","CTSD","CTSL",
            "ATP6V0D1","ATP6V1A")
cargo  <- c("INS","CHGA","CHGB","PCSK1","PCSK2","CPE",
            "RAB3A","RAB27A","CADPS","SYTL4","PCLO","RAB3GAP1")

cat("\n===== DESCRIPTIVE: fusion vs cargo (NO TEST) =====\n")
for (mm in list(list(m1,"M1"), list(m2,"M2"))) {
  r <- mm[[1]]
  cat(sprintf("  %s  fusion median = %.3f (n=%d)   cargo median = %.3f (n=%d)\n",
      mm[[2]],
      median(r$sd_ratio[r$gene %in% fusion], na.rm=TRUE),
      sum(r$gene %in% fusion),
      median(r$sd_ratio[r$gene %in% cargo], na.rm=TRUE),
      sum(r$gene %in% cargo)))
}
cat("\n  This grouping was suggested BY the results and is not part of\n")
cat("  the frozen gene set. It is stated as a hypothesis for\n")
cat("  independent testing, not as a finding. Do not compute a\n")
cat("  p-value for it in this dataset.\n")


# ---- 8. Figure ------------------------------------------------------

lab <- c("UNC13D","STX7","GABARAP","INS","STX1A","RAB3A","RAB27A","CTSB")

p <- m2 %>%
  mutate(highlight = gene %in% lab) %>%
  ggplot(aes(x = sd_ratio, y = reorder(gene, sd_ratio), colour = set)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = quantile(bg2, .99), linetype = "dotted",
             colour = "grey30") +
  geom_pointrange(aes(xmin = ci_lo, xmax = ci_hi), size = .3) +
  scale_colour_manual(values = c(crinophagy = "#1b5e20",
                                 comparator = "#90a4ae",
                                 `post-hoc` = "#ef6c00")) +
  labs(x = "SD ratio T2D / ND (95% bootstrap CI)", y = NULL, colour = NULL,
       title = "Variance ratios, fully adjusted model",
       subtitle = "Dotted line: 99th percentile of 20,255-gene background") +
  theme_minimal(base_size = 9)

ggsave(file.path(REP, "fig_canonical_forest.pdf"), p, width = 7, height = 9)

cat("\nWrote canonical_both_models.csv, canonical_full.csv,",
    "fig_canonical_forest.pdf\n")
cat("\nFREEZE HERE. No further gene panels, filters, or backgrounds.\n")
