# =====================================================================
# 55_open_clinical.R
#
# Resolves the computable parts of four OPEN audit items:
#   R2-M6  21-donor diagnostic reclassification sensitivity;
#          retained replication donor table and batch x status
#   R2-M3  mean effects with CIs in the same matrices as the variance
#          analysis; insulin-content estimates
#   R2-m4  the CI 0.97-2.21 against Brown-Forsythe P = 0.018
#   R2-M10 core-fusion-machinery sensitivity  (STOPS until you classify
#          the 31 genes; it is not my call to make)
#
# RUN:  source("40_common.R"); source("55_open_clinical.R")
# =====================================================================

stopifnot(exists("mat"), exists("pd2"), exists("crino"), exists("macro"))
set.seed(1)
grp <- factor(as.character(pd2$group), levels = c("ND", "T2D"))

TARG <- c(UNC13D = "ENSG00000092929", STX7 = "ENSG00000079950",
          GABARAP = "ENSG00000170296", CTSD = "ENSG00000117984",
          RAB3A = "ENSG00000105649")
TARG <- TARG[TARG %in% rownames(mat)]


# =====================================================================
# 1 | R2-M6a  the 21 reclassified discovery donors
# =====================================================================

cat("\n############ 1. DIAGNOSTIC RECLASSIFICATION (R2-M6) ############\n")
don <- read.csv("donor.csv", stringsAsFactors = FALSE)
P   <- readRDS("proteomics/proc_prot.rds")
Mp  <- as.matrix(P[, setdiff(names(P), "gene_id")]); mode(Mp) <- "numeric"
rownames(Mp) <- as.character(P$gene_id)
dd  <- don[match(colnames(Mp), don$record_id), ]

diff_ix <- which(dd$diagnosis != dd$diagnosis_computed)
cat(sprintf("donors whose chart diagnosis differs from diagnosis_computed: %d\n",
            length(diff_ix)))
if (length(diff_ix))
  print(table(chart = dd$diagnosis[diff_ix],
              computed = dd$diagnosis_computed[diff_ix]))

COV <- c("donorage","donorsex","bodymassindex",
         "predistributionculturetime","percentieqrecoverypostculture")
DISC_ENTREZ <- c(UNC13D = "201294", RAB3A = "5864", GABARAP = "11337",
                 CTSD = "1509")

disc_sd <- function(label_col) {
  g <- ifelse(dd[[label_col]] == "Type2", "T2D",
       ifelse(dd[[label_col]] == "Pre.T2D", "preT2D", "ND"))
  ok <- complete.cases(dd[, COV]) & g %in% c("ND","T2D")
  d  <- dd[ok, ]; d$donorsex <- factor(d$donorsex); gg <- g[ok]
  X  <- model.matrix(~ donorage + donorsex + bodymassindex +
                       predistributionculturetime +
                       percentieqrecoverypostculture, data = d)
  Y  <- t(Mp[, ok, drop = FALSE])
  R  <- lm.fit(X, Y)$residuals
  cs <- function(M) sqrt(colSums(sweep(M, 2, colMeans(M))^2) / (nrow(M) - 1))
  r  <- cs(R[gg == "T2D", , drop = FALSE]) / cs(R[gg == "ND", , drop = FALSE])
  list(ratios = r, n = table(gg))
}

a <- disc_sd("diagnosis_computed")
b <- disc_sd("diagnosis")
cat("\ngroup sizes, computed diagnosis:\n"); print(a$n)
cat("group sizes, chart diagnosis:\n");    print(b$n)

rc <- do.call(rbind, lapply(names(DISC_ENTREZ), function(s) {
  e <- DISC_ENTREZ[[s]]
  data.frame(gene = s,
             computed_diagnosis = if (e %in% names(a$ratios)) a$ratios[[e]] else NA,
             chart_diagnosis    = if (e %in% names(b$ratios)) b$ratios[[e]] else NA)
}))
rc$change <- rc$chart_diagnosis - rc$computed_diagnosis
cat("\ndiscovery SD ratios under each diagnostic definition:\n")
print(rc, row.names = FALSE, digits = 4)
cat("\n  Report this table. A result that moves materially between the two\n")
cat("  definitions is a result that depends on a classification choice.\n")
write.csv(rc, "C_reclassification_sensitivity.csv", row.names = FALSE)


# =====================================================================
# 2 | R2-M6b  retained replication donors: batch x status, clinical fields
# =====================================================================

cat("\n############ 2. RETAINED REPLICATION DONORS (R2-M6) ############\n")
cat("batch x glycemic status, the 54 analysed donors:\n")
print(table(batch = pd2$batch, status = pd2$group))
cat("\nchi-square test of association:\n")
print(suppressWarnings(chisq.test(table(pd2$batch, pd2$group))))

meta_f <- "replication_GSE164416/gse164416_metadata.csv"
if (file.exists(meta_f)) {
  md <- read.csv(meta_f, stringsAsFactors = FALSE)
  cat(sprintf("\n%s: %d rows x %d cols\n", meta_f, nrow(md), ncol(md)))
  print(names(md))
  path_cols <- grep("patho|tumou|tumor|malign|diagn|surg|indicat|histol",
                    names(md), ignore.case = TRUE, value = TRUE)
  if (length(path_cols)) {
    cat("\npathology-related columns found:\n")
    for (p in path_cols) { cat("\n--", p, "\n"); print(table(md[[p]], useNA = "ifany")) }
  } else {
    cat("\nNo pathology column in the GEO metadata. The tumour counts the\n")
    cat("reviewer cites (27/39 and 9/18) are from the SOURCE cohort, not\n")
    cat("these 54 donors. Say that explicitly, and if the retained-donor\n")
    cat("pathology is not available, say that too rather than implying it is.\n")
  }
} else cat("\nMISSING:", meta_f, "\n")


# =====================================================================
# 3 | R2-M3  mean effects with CIs, and insulin content
# =====================================================================

cat("\n############ 3. MEAN EFFECTS WITH CIs (R2-M3) ############\n")
me <- do.call(rbind, lapply(names(TARG), function(s) {
  e <- TARG[[s]]
  d <- pd2; d$y <- mat[e, rownames(pd2)]
  f <- lm(y ~ group + batch + imm + beta, data = d)
  ci <- confint(f)["groupT2D", ]
  data.frame(gene = s, mean_diff = coef(f)[["groupT2D"]],
             lo = ci[1], hi = ci[2],
             p = summary(f)$coefficients["groupT2D", 4],
             sd_ratio = unname(sd_ratios(RES2, as.character(grp))[e]))
}))
cat("replication cohort, same matrix and model as the variance analysis:\n")
print(me, row.names = FALSE, digits = 4)
cat("\n  Quote the interval, not the P value. 'No difference detected' with\n")
cat("  an interval spanning a meaningful effect is not evidence of\n")
cat("  preserved abundance, which is exactly R2-M3's point.\n")
write.csv(me, "C_mean_effects_with_CI.csv", row.names = FALSE)

cat("\n-- insulin content (discovery donors) --\n")
ic <- c("pdinsulincontent","pdinsulinperieq","pdinsulindnaratio")
ic <- ic[ic %in% names(dd)]
for (v in ic) {
  g <- ifelse(dd$diagnosis_computed == "Type2", "T2D", "non-T2D")
  x <- dd[[v]]
  if (all(is.na(x))) { cat(sprintf("  %-20s all missing\n", v)); next }
  tt <- try(t.test(x ~ g), silent = TRUE)
  if (inherits(tt, "try-error")) next
  cat(sprintf("  %-20s non-T2D %.3g vs T2D %.3g   diff CI [%.3g, %.3g]  P = %.3g\n",
              v, tt$estimate[1], tt$estimate[2],
              tt$conf.int[1], tt$conf.int[2], tt$p.value))
}


# =====================================================================
# 4 | R2-m4  the confidence interval against Brown-Forsythe
# =====================================================================

cat("\n############ 4. CI vs BROWN-FORSYTHE (R2-m4) ############\n")
e <- TARG[["UNC13D"]]
r <- RES2[, e]; a <- r[grp == "T2D"]; b <- r[grp == "ND"]
bt <- replicate(20000, sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))
q  <- quantile(bt, c(.025, .975))
z  <- abs(r - ave(r, grp, FUN = median))
bf <- t.test(z[grp == "T2D"], z[grp == "ND"], var.equal = TRUE)

cat(sprintf("  SD ratio          : %.3f   bootstrap 95%% CI [%.3f, %.3f]\n",
            sd(a)/sd(b), q[1], q[2]))
cat(sprintf("  Brown-Forsythe    : t = %.3f, P = %.4f\n", bf$statistic, bf$p.value))
cat(sprintf("  mean |deviation|  : T2D %.3f vs ND %.3f, difference CI [%.3f, %.3f]\n",
            mean(z[grp=="T2D"]), mean(z[grp=="ND"]),
            -bf$conf.int[2], -bf$conf.int[1]))
cat(sprintf("  bootstrap skewness: %.2f  (median %.3f vs point estimate %.3f)\n",
            mean((bt - mean(bt))^3) / sd(bt)^3, median(bt), sd(a)/sd(b)))
cat("\n  These are different estimands, not two readings of one quantity:\n")
cat("  the CI is for a RATIO of SDs, right-skewed at this sample size; BF\n")
cat("  tests a DIFFERENCE in mean absolute deviation. Report both and say\n")
cat("  so. Do not present the significant one as the answer.\n")


# =====================================================================
# 5 | R2-M10  core fusion machinery  -- NEEDS YOUR CLASSIFICATION
# =====================================================================

cat("\n############ 5. CORE-FUSION SENSITIVITY (R2-M10) ############\n")
gs <- read.delim("frozen_gene_sets_v3.txt", stringsAsFactors = FALSE)
cr <- gs[gs$set == "crinophagy", ]
cat("the 31 prospectively defined crinophagy genes:\n")
print(cr[, intersect(names(cr), c("SYMBOL","gene","ENTREZID","set"))],
      row.names = FALSE)

CORE_SYMBOLS <- c("UNC13D","STX7","STX8","VAMP7","VAMP8","SNAP23",
                  "STX17","VPS11","VPS16","VPS33A","RAB27A","SYTL4")

if (!length(CORE_SYMBOLS)) {
  cat("\n*** CORE_SYMBOLS is empty; classify the 31 genes and source again.\n")
} else {
  ## frozen_gene_sets_v3.txt is keyed by ENTREZID, not ENSEMBL, so map first
  ## (1:many resolved by highest mean expression, as in scripts 33 and 35)
  suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
  mm <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(cr$ENTREZID),
                              keytype = "ENTREZID", columns = "ENSEMBL")
  detach("package:org.Hs.eg.db", unload = TRUE)
  detach("package:AnnotationDbi", unload = TRUE)
  mm <- mm[!is.na(mm$ENSEMBL) & mm$ENSEMBL %in% names(bg2), ]
  mm$expr <- rowMeans(mat[mm$ENSEMBL, , drop = FALSE])
  mm <- mm[order(mm$ENTREZID, -mm$expr), ]
  mm <- mm[!duplicated(mm$ENTREZID), ]
  cr$ENSEMBL <- mm$ENSEMBL[match(as.character(cr$ENTREZID), mm$ENTREZID)]

  sym2ens <- setNames(cr$ENSEMBL, cr$SYMBOL)
  missing <- CORE_SYMBOLS[is.na(sym2ens[CORE_SYMBOLS])]
  if (length(missing))
    cat("not mapped to the replication matrix:", paste(missing, collapse = ", "), "\n")

  set_test <- function(symbols, label, nperm = 1e5) {
    g <- intersect(unname(sym2ens[symbols]), names(bg2))
    g <- g[!is.na(g)]
    if (!length(g)) { cat("  ", label, ": no genes mapped\n"); return(invisible()) }
    obs <- median(bg2[g])
    nul <- replicate(nperm, median(sample(bg2, length(g))))
    cat(sprintf("  %-34s n=%2d  median %.4f  null %.4f  P = %.5f\n",
                label, length(g), obs, median(nul),
                (sum(nul >= obs) + 1) / (nperm + 1)))
  }

  cat("\n-- set-level comparison, Model 2 (fully adjusted) --\n")
  set_test(cr$SYMBOL,                              "Full frozen set (31)")
  set_test(CORE_SYMBOLS,                           "Core fusion machinery")
  set_test(setdiff(CORE_SYMBOLS, "STX17"),         "Core without STX17")
  set_test(setdiff(cr$SYMBOL, CORE_SYMBOLS),       "Non-core remainder")

  cat("\n  Report alongside the full set, never in place of it. The frozen\n")
  cat("  set is unchanged; this is a post-hoc functional sensitivity analysis\n")
  cat("  and must be labelled as one. STX17 is shown separately because its\n")
  cat("  membership in a SG-lysosome fusion core is contestable.\n")
}

cat("\n#################### DONE ####################\n")
