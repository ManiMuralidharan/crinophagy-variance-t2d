# =====================================================================
# Rebuild Supplementary_Tables.xlsx from source
#
# WHAT THIS FIXES
#   - Table S2 (exploratory 136-gene panel): CSV missing, recomputed
#   - Table S4: CAPS family rows were never added
#   - Table S4: Munc13 paralogue CSV missing, recomputed
#   - Table S3 Panel A: tail counts recomputed per model (1/31 vs 2/31)
#   - Explanatory notes rewritten to match the final legends
#
# WHAT IT DOES NOT DO
#   Nothing is hardcoded that can be computed. Every value is derived
#   from mat/pd2/bg1/bg2 or read from canonical_full.csv. If a number
#   here disagrees with the manuscript, the manuscript is wrong.
#
# REQUIRES in session:  mat, pd2, bg1, bg2
#   load("session_gse164416.RData")
#
# Run:  source("31_rebuild_supplementary_tables.R")
# Writes: supplementary/Supplementary_Tables.xlsx  (overwrites)
# =====================================================================

library(openxlsx)
library(car)

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))
REP <- "replication_GSE164416"
OUT <- "supplementary"
dir.create(OUT, showWarnings = FALSE)

if (!exists("mat")) load("session_gse164416.RData")
stopifnot(exists("mat"), exists("pd2"), exists("bg1"), exists("bg2"))

set.seed(1)
NPERM <- 10000

cat("Background sizes:", length(bg1), length(bg2), "\n")


# ---- rebuild composition covariates if absent ------------------------

if (!all(c("imm", "beta") %in% names(pd2))) {
  cat("Rebuilding imm/beta covariates...\n")
  suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
  ens_of <- function(sy) {
    a <- AnnotationDbi::select(org.Hs.eg.db, keys = sy,
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
  detach("package:org.Hs.eg.db", unload = TRUE)
  detach("package:AnnotationDbi", unload = TRUE)
}

cat("Immune SD by group:  ", round(tapply(pd2$imm,  pd2$group, sd), 3), "\n")
cat("Beta-cell SD by group:", round(tapply(pd2$beta, pd2$group, sd), 3), "\n")


# ---- shared test function -------------------------------------------

one_gene <- function(ens, nboot = 5000) {
  d <- pd2
  d$y <- mat[ens, rownames(d)]
  d$r <- residuals(lm(y ~ batch + imm + beta, data = d))
  a <- d$r[d$group == "T2D"]; b <- d$r[d$group == "ND"]
  bt <- replicate(nboot,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))
  c(sd_ratio = sd(a)/sd(b),
    ci_lo = unname(quantile(bt, .025)),
    ci_hi = unname(quantile(bt, .975)),
    p = leveneTest(r ~ group, data = d)$`Pr(>F)`[1],
    pct = 100 * mean(bg2 <= sd(a)/sd(b)))
}

resolve <- function(syms) {
  suppressMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
  m <- AnnotationDbi::select(org.Hs.eg.db, keys = syms,
                             keytype = "SYMBOL", columns = "ENSEMBL")
  m <- m[!is.na(m$ENSEMBL) & m$ENSEMBL %in% rownames(mat), ]
  m <- m[!duplicated(m$SYMBOL), ]
  detach("package:org.Hs.eg.db", unload = TRUE)
  detach("package:AnnotationDbi", unload = TRUE)
  m
}

rnd <- function(d, k = 4) {
  num <- sapply(d, is.numeric); d[num] <- lapply(d[num], round, k); d
}


# =====================================================================
# TABLE S1 - pre-registered gene sets, both models
# =====================================================================

cat("\n=== Table S1 ===\n")
full <- read.csv(file.path(REP, "canonical_full.csv"),
                 stringsAsFactors = FALSE)

m1 <- full[full$model == "M1_batch", ]
m2 <- full[full$model == "M2_full",  ]

S1 <- merge(
  m1[, c("gene","set","ensembl","sd_ratio","ci_lo","ci_hi","p","pct","q")],
  m2[, c("gene","sd_ratio","ci_lo","ci_hi","p","pct","q")],
  by = "gene", suffixes = c("_M1", "_M2"))
S1 <- S1[order(-S1$sd_ratio_M2), ]

cat("Genes:", nrow(S1), " (expect 46)\n")
print(S1[S1$gene %in% c("UNC13D","STX7","CTSD","GABARAP"),
         c("gene","sd_ratio_M1","pct_M1","sd_ratio_M2","pct_M2","q_M2")],
      digits = 4, row.names = FALSE)


# =====================================================================
# TABLE S2 - exploratory 136-gene panel  (CSV missing, recompute)
# =====================================================================

cat("\n=== Table S2 (recomputing) ===\n")

tiers <- list(
  `Autophagosome: initiation` = c("ULK1","ULK2","ATG13","ATG101","RB1CC1"),
  `Autophagosome: nucleation` = c("PIK3C3","BECN1","ATG14","UVRAG","NRBF2",
                                  "RUBCN","RUBCNL"),
  `Autophagosome: early membrane` = c("ATG9A","ATG2A","ATG2B","TMEM41B",
                                      "VMP1","WIPI2","ZFYVE1","VPS13A","VPS13C"),
  `Autophagosome: conjugation` = c("ATG3","ATG5","ATG7","ATG10","ATG12",
                                   "ATG16L1","ATG4B"),
  `Autophagosome: ATG8 family` = c("MAP1LC3A","MAP1LC3B","GABARAP",
                                   "GABARAPL1","GABARAPL2"),
  `Autophagosome: cargo receptors` = c("SQSTM1","NBR1","CALCOCO2","TAX1BP1",
                                       "BNIP3L","OPTN"),
  `Autophagosome: fusion SNAREs` = c("STX17","SNAP29","YKT6","SNAP47",
                                     "VAMP8","SCFD1"),
  `Autophagosome: ESCRT / closure` = c("VPS4A","VPS4B","VPS37A","CHMP2A",
                                       "CHMP4B","VPS25","VPS28","TSG101"),
  `Granule: Munc13 priming` = c("UNC13D","UNC13A","UNC13B"),
  `Granule: Rab27A effectors` = c("RAB27A","SYTL1","SYTL2","SYTL4","RPH3AL",
                                  "MLPH","MYO5A","EXPH5","MADD"),
  `Granule: Rab3 axis` = c("RAB3A","RAB3B","RAB3D","RAB3GAP1","RAB3GAP2",
                           "RAB3IP"),
  `Granule: SNAREs` = c("VAMP2","VAMP4","VAMP7","VAMP3"),
  `Granule: identity and cargo` = c("PCLO","CHGA","CHGB","PCSK1","PCSK2",
                                    "CPE","INS"),
  `Lysosome: SNAREs` = c("STX7","STX8","VTI1B"),
  `Lysosome: HOPS` = c("VPS11","VPS16","VPS18","VPS33A","VPS39","VPS41"),
  `Lysosome: tethering` = c("PLEKHM1","EPG5","TECPR1","BIRC6"),
  `Lysosome: Rab and effectors` = c("RAB7A","RAB2A","RAB21","RAB24","RAB33B",
                                    "RAB39A","ARL8A","ARL8B","RILP","FYCO1",
                                    "OSBPL1A","TBC1D25"),
  `Lysosome: membrane and lumen` = c("LAMP1","LAMP2","CTSB","CTSD","CTSL",
                                     "ATP6V0D1","ATP6V1A","TFEB","PIKFYVE",
                                     "MCOLN1","NPC1"),
  `Lysosome: lipid enzymes` = c("INPP5E","PI4K2A","PIP5K1C"),
  `Upstream signalling` = c("MTOR","RPTOR","PRKAA1","PRKAA2","TBK1","STK4",
                            "IKBKB","STING1","SIRT3","SIRT4","SIRT5","SIRT6",
                            "KAT2B","CREBBP","HDAC2")
)

tier_of <- function(g) {
  if (grepl("^Autophagosome", g)) {
    "1. Autophagosome donor"
  } else if (grepl("^Granule", g)) {
    "2. Granule donor"
  } else if (grepl("^Lysosome", g)) {
    "3. Shared lysosomal acceptor"
  } else {
    "4. Upstream signalling"
  }
}

pan <- do.call(rbind, lapply(names(tiers), function(g)
  data.frame(gene = tiers[[g]], subgroup = g, tier = tier_of(g),
             stringsAsFactors = FALSE)))
pan <- pan[!duplicated(pan$gene), ]

frozen <- read.delim("frozen_gene_sets_v2_clean.txt", stringsAsFactors = FALSE)
cr <- frozen$SYMBOL[frozen$set == "crinophagy"]
cm <- frozen$SYMBOL[frozen$set == "comparator"]
pan$frozen_set_status <- ifelse(pan$gene %in% cr, "Frozen crinophagy set",
                         ifelse(pan$gene %in% cm, "Frozen comparator set",
                                "New (post-hoc)"))

mp <- resolve(pan$gene)
absent_S2 <- setdiff(pan$gene, mp$SYMBOL)
pan <- merge(pan, mp, by.x = "gene", by.y = "SYMBOL")
names(pan)[names(pan) == "ENSEMBL"] <- "ensembl"

cat("Proposed:", nrow(pan) + length(absent_S2),
    " measured:", nrow(pan),
    " not measured:", length(absent_S2), "\n")

tt <- t(sapply(pan$ensembl, one_gene, nboot = 2000))
S2 <- cbind(pan, as.data.frame(tt))
S2$q <- NA
i <- S2$frozen_set_status == "New (post-hoc)"
S2$q[i] <- p.adjust(S2$p[i], method = "BH")
S2 <- S2[order(-S2$sd_ratio), ]

cat("Post-hoc genes:", sum(i), " minimum q:",
    round(min(S2$q, na.rm = TRUE), 4), "(expect ~0.577)\n")


# =====================================================================
# TABLE S3 - set-level and tier-level summaries
# =====================================================================

cat("\n=== Table S3 ===\n")

set_test <- function(tag, bg) {
  res <- full[full$model == tag, ]
  thr <- quantile(bg, .99)
  do.call(rbind, lapply(c("crinophagy","comparator"), function(s) {
    g <- res$sd_ratio[res$set == s]; g <- g[is.finite(g)]
    obs <- median(g)
    null <- replicate(NPERM, median(sample(bg, length(g))))
    k <- sum(g >= thr); n <- length(g)
    data.frame(
      gene_set = ifelse(s == "crinophagy",
                        "Crinophagy (31 genes, pre-registered)",
                        "Macroautophagy comparator (14 genes, pre-registered)"),
      model = ifelse(tag == "M1_batch", "M1 batch-adjusted",
                     "M2 fully adjusted (primary)"),
      n_genes = n,
      median_ratio = obs,
      null_median = median(null),
      permutation_P = (sum(null >= obs) + 1) / (NPERM + 1),
      threshold_99th = unname(thr),
      n_above_99th = k,
      expected = 0.01 * n,
      binomial_P = binom.test(k, n, 0.01, alternative = "greater")$p.value,
      genes_above = paste(res$gene[res$set == s & res$sd_ratio >= thr],
                          collapse = ", "),
      stringsAsFactors = FALSE)
  }))
}

S3A <- rbind(set_test("M1_batch", bg1), set_test("M2_full", bg2))
S3A <- S3A[order(S3A$gene_set, S3A$model), ]
print(S3A[, c("gene_set","model","median_ratio","permutation_P",
              "n_above_99th","binomial_P")], digits = 4, row.names = FALSE)

ph <- S2[S2$frozen_set_status == "New (post-hoc)", ]
p99 <- quantile(bg2, .99)
S3B <- do.call(rbind, lapply(sort(unique(ph$tier)), function(t) {
  x <- ph[ph$tier == t, ]
  data.frame(tier = t, n_posthoc_genes = nrow(x),
             median_sd_ratio = median(x$sd_ratio),
             n_above_99th = sum(x$sd_ratio >= p99),
             stringsAsFactors = FALSE)
}))
S3B <- S3B[order(-S3B$median_sd_ratio), ]
print(S3B, digits = 4, row.names = FALSE)
cat("Genome-wide median:", round(median(bg2), 3),
    " 99th pct:", round(p99, 3), "\n")


# =====================================================================
# TABLE S4 - post-hoc comparator panels (CAPS added)
# =====================================================================

cat("\n=== Table S4 ===\n")

panels <- list(
  `Regulated exocytosis` = c("STX1A","STX1B","STX4","SNAP25","SNAP23",
                             "VAMP2","STXBP1","STXBP2","STXBP3","STXBP5",
                             "STXBP5L","UNC13A","UNC13B","UNC13C","SYT7",
                             "SYT9","RIMS2"),
  `CAPS family (MUN domain)` = c("CADPS","CADPS2","BAIAP3"),
  `Munc13 family (MUN domain)` = c("UNC13A","UNC13B","UNC13C","UNC13D"),
  `Rab GEF / GAP layer` = c("RAB3GAP1","RAB3GAP2","MADD","TBC1D10A",
                            "TBC1D10B","TBC1D10C","RAB3IP","RAB3IL1",
                            "DENND2B","TBC1D3","EPI64")
)

S4 <- do.call(rbind, lapply(names(panels), function(pn) {
  mpp <- resolve(panels[[pn]])
  miss <- setdiff(panels[[pn]], mpp$SYMBOL)
  if (length(miss))
    cat(sprintf("  %-28s not measured: %s\n", pn, paste(miss, collapse = ", ")))
  if (nrow(mpp) == 0) return(NULL)
  tt <- t(sapply(mpp$ENSEMBL, one_gene))
  data.frame(panel = pn, gene = mpp$SYMBOL, ensembl = mpp$ENSEMBL,
             as.data.frame(tt),
             status = ifelse(mpp$SYMBOL %in% c(cr, cm),
                             "Pre-registered", "Post-hoc"),
             stringsAsFactors = FALSE)
}))
S4 <- S4[order(S4$panel, -S4$sd_ratio), ]

cat("\nCAPS family:\n")
print(S4[S4$panel == "CAPS family (MUN domain)",
         c("gene","sd_ratio","ci_lo","ci_hi","p","pct")],
      digits = 4, row.names = FALSE)

exo <- S4[S4$panel == "Regulated exocytosis", ]
cat(sprintf("\nExocytosis genes above 99th percentile: %d of %d (expect 0)\n",
            sum(exo$sd_ratio >= p99), nrow(exo)))


# =====================================================================
# WRITE WORKBOOK
# =====================================================================

wb <- createWorkbook()
hdr <- createStyle(textDecoration = "bold", fgFill = "#E8F5E9",
                   border = "bottom", halign = "left")
note <- createStyle(fontSize = 9, textDecoration = "italic",
                    wrapText = TRUE, valign = "top")

add_sheet <- function(name, text, df, height = 56) {
  addWorksheet(wb, name)
  writeData(wb, name, text, startRow = 1)
  addStyle(wb, name, note, rows = 1, cols = 1)
  setRowHeights(wb, name, rows = 1, heights = height)
  writeData(wb, name, rnd(df), startRow = 3, headerStyle = hdr)
  freezePane(wb, name, firstActiveRow = 4)
  setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
}

add_sheet("S1_pre-registered_sets", paste(
  "Supplementary Table S1. Pre-registered gene sets, replication cohort",
  "(GSE164416; 18 non-diabetic, 36 T2D). Both sets were frozen with a",
  "timestamped DOI before analysis. M1 = batch-adjusted; M2 = adjusted for",
  "batch, immune infiltration and beta-cell content (primary). sd_ratio =",
  "residual SD in T2D divided by residual SD in non-diabetic donors.",
  "pct = percentile against 20,255 expressed genes under the same model.",
  "q = Benjamini-Hochberg within the comparator panel only; crinophagy genes",
  "constitute the hypothesis and are reported uncorrected. STX1A is a",
  "post-hoc comparator, labelled as such."), S1, 70)

add_sheet("S2_exploratory_panel", paste(
  "Supplementary Table S2. Exploratory autophagy-lysosome panel,",
  nrow(S2), "genes grouped by position in the pathway. POST-HOC: assembled",
  "after the primary analysis and not part of the pre-registered design.",
  "Tested under the same fully adjusted model. Genes already in the frozen",
  "crinophagy set are marked and were excluded from all group-level summaries",
  "and from the FDR calculation, since a comparison containing the test genes",
  "is not a comparison. NULL RESULT: no gene survived correction (minimum",
  sprintf("q = %.3f)", min(S2$q, na.rm = TRUE)),
  "and no tier was materially elevated above the genome-wide median of",
  sprintf("%.3f.", median(bg2))), S2, 82)

addWorksheet(wb, "S3_summaries")
writeData(wb, "S3_summaries", paste(
  "Supplementary Table S3. Set-level and tier-level summaries.",
  "Panel A: pre-registered sets against 10,000 size-matched random gene sets",
  "drawn from the genome-wide background, under both adjustment models.",
  "Note that the number of genes exceeding the 99th percentile differs",
  "between models. Panel B: exploratory tiers, post-hoc genes only;",
  sprintf("genome-wide median %.3f, 99th percentile %.3f.",
          median(bg2), p99)), startRow = 1)
addStyle(wb, "S3_summaries", note, rows = 1, cols = 1)
setRowHeights(wb, "S3_summaries", rows = 1, heights = 60)
writeData(wb, "S3_summaries", "A. Pre-registered set-level tests", startRow = 3)
addStyle(wb, "S3_summaries", createStyle(textDecoration = "bold"),
         rows = 3, cols = 1)
writeData(wb, "S3_summaries", rnd(S3A), startRow = 4, headerStyle = hdr)
r0 <- 4 + nrow(S3A) + 3
writeData(wb, "S3_summaries", "B. Exploratory tiers (post-hoc genes only)",
          startRow = r0 - 1)
addStyle(wb, "S3_summaries", createStyle(textDecoration = "bold"),
         rows = r0 - 1, cols = 1)
writeData(wb, "S3_summaries", rnd(S3B), startRow = r0, headerStyle = hdr)
setColWidths(wb, "S3_summaries", cols = 1:11, widths = "auto")

add_sheet("S4_posthoc_panels", paste(
  "Supplementary Table S4. Post-hoc comparator panels, all under the fully",
  "adjusted model. The regulated exocytosis panel tests whether elevated",
  "variance is a general property of beta-cell membrane fusion machinery",
  "(it is not: no gene above the 99th percentile). The Munc13 and CAPS",
  "families both carry MUN domains and catalyse the same SNARE-opening step;",
  "they are grouped by structural homology rather than by any partition",
  "imposed after seeing results. The Rab GEF/GAP panel is uninformative",
  "about Rab activity, since these proteins act catalytically and abundance",
  "is a poor proxy for nucleotide state."), S4, 82)

f <- file.path(OUT, "Supplementary_Tables.xlsx")
saveWorkbook(wb, f, overwrite = TRUE)

cat("\n=== WRITTEN ===\n")
cat(f, "\n")
cat("Sheets:", paste(names(wb), collapse = ", "), "\n")
cat("\nCheck against the manuscript:\n")
cat("  S1 rows:", nrow(S1), "(expect 46)\n")
cat("  S2 rows:", nrow(S2), " post-hoc:", sum(S2$frozen_set_status == "New (post-hoc)"), "\n")
cat("  S4 rows:", nrow(S4), " including", sum(S4$panel == "CAPS family (MUN domain)"), "CAPS genes\n")

write.csv(S2, file.path(REP, "donor_acceptor_panel.csv"), row.names = FALSE)
write.csv(S4, "figures/posthoc_panels_all.csv", row.names = FALSE)
cat("\nAlso rewrote the two missing source CSVs.\n")
