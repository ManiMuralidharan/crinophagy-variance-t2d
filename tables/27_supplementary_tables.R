# =====================================================================
# Supplementary tables
#
#   S1  Pre-registered gene sets, both models, full results
#   S2  Exploratory autophagy-lysosome panel (136 genes, null result)
#   S3  Set-level and tier-level summaries
#   S4  Post-hoc comparator panels (exocytosis, Rab regulators, Munc13)
#
# Written as one .xlsx workbook so a reader can check every number
# behind the manuscript, including the negative results.
#
# Requires: res, tab, ph in session (from the extended panel), plus
# canonical_full.csv and the panel CSVs on disk.
# =====================================================================

setwd("/Users/mmani3/Downloads/Human Islets")
REP  <- "replication_GSE164416"
FIGS <- "figures"
OUT  <- "supplementary"
dir.create(OUT, showWarnings = FALSE)

if (!requireNamespace("openxlsx", quietly = TRUE))
  install.packages("openxlsx")
library(openxlsx)

rnd <- function(d, digits = 4) {
  num <- sapply(d, is.numeric)
  d[num] <- lapply(d[num], round, digits)
  d
}

wb <- createWorkbook()

hdr <- createStyle(textDecoration = "bold", fgFill = "#E8F5E9",
                   border = "bottom", halign = "left")

add_sheet <- function(name, note, df) {
  addWorksheet(wb, name)
  writeData(wb, name, note, startRow = 1)
  addStyle(wb, name, createStyle(fontSize = 9, textDecoration = "italic",
                                 wrapText = TRUE),
           rows = 1, cols = 1)
  setRowHeights(wb, name, rows = 1, heights = 46)
  writeData(wb, name, df, startRow = 3, headerStyle = hdr)
  freezePane(wb, name, firstActiveRow = 4)
  setColWidths(wb, name, cols = 1:ncol(df), widths = "auto")
}


# ---- S1: pre-registered sets ----------------------------------------

full <- read.csv(file.path(REP, "canonical_full.csv"),
                 stringsAsFactors = FALSE)

s1 <- reshape(
  full[, c("gene","set","model","sd_ratio","ci_lo","ci_hi","p","pct","q")],
  idvar = c("gene","set"), timevar = "model", direction = "wide")

names(s1) <- sub("\\.M1_batch$", "_M1", names(s1))
names(s1) <- sub("\\.M2_full$",  "_M2", names(s1))
s1 <- s1[order(-s1$sd_ratio_M2), ]

add_sheet("S1_frozen_sets",
  paste("Supplementary Table S1. Pre-registered gene sets, replication cohort",
        "(GSE164416; 18 non-diabetic, 36 T2D). Both sets were frozen with a",
        "timestamped DOI before analysis. M1 = batch-adjusted; M2 = adjusted for",
        "batch, immune infiltration and beta-cell content (primary). sd_ratio =",
        "residual SD in T2D divided by residual SD in non-diabetic donors.",
        "pct = percentile against 20,255 expressed genes under the same model.",
        "q = Benjamini-Hochberg within the comparator panel only; crinophagy genes",
        "constitute the hypothesis and are reported uncorrected. STX1A is a post-hoc",
        "comparator, labelled as such."),
  rnd(s1))


# ---- S2: exploratory panel (the null result) -------------------------

if (exists("res")) {
  s2 <- res[, c("SYMBOL","tier","subgroup","status","ratio","lo","hi",
                "p","pct","q")]
  names(s2) <- c("gene","tier","subgroup","frozen_set_status",
                 "sd_ratio","ci_lo","ci_hi","p","percentile","q_BH")
  s2 <- s2[order(-s2$sd_ratio), ]
} else if (file.exists(file.path(FIGS, "donor_acceptor_panel.csv"))) {
  s2 <- read.csv(file.path(FIGS, "donor_acceptor_panel.csv"),
                 stringsAsFactors = FALSE)
} else stop("Extended panel result not found.")

add_sheet("S2_exploratory_panel",
  paste("Supplementary Table S2. Exploratory autophagy-lysosome panel, 136 genes",
        "grouped by position in the pathway. POST-HOC: assembled after the primary",
        "analysis and not part of the pre-registered design. Tested under the same",
        "fully adjusted model. Genes already in the frozen crinophagy set are marked",
        "and were excluded from all group-level summaries and from the FDR",
        "calculation, since a comparison containing the test genes is not a",
        "comparison. NULL RESULT: no gene survived correction (minimum q = 0.577)",
        "and no tier was materially elevated above the genome-wide median of 1.077."),
  rnd(s2))


# ---- S3: summaries ---------------------------------------------------

set_level <- data.frame(
  analysis = c("Crinophagy set (31 genes)", "Crinophagy set (31 genes)",
               "Comparator set (14 genes)", "Comparator set (14 genes)"),
  model = c("M1 batch-adjusted", "M2 fully adjusted",
            "M1 batch-adjusted", "M2 fully adjusted"),
  median_ratio = c(1.228, 1.178, 1.192, 1.036),
  null_median  = c(1.089, 1.079, 1.091, 1.080),
  permutation_p = c(0.0066, 0.036, 0.103, 0.719),
  n_above_99th = c(1, 2, 1, 1),
  expected     = c(0.31, 0.31, 0.14, 0.14),
  binomial_p   = c(0.268, 0.038, 0.131, 0.131),
  stringsAsFactors = FALSE)

if (exists("tab")) {
  tier_level <- tab[order(-tab$ratio), c("tier","n","ratio","above99")]
  names(tier_level) <- c("tier","n_posthoc_genes","median_sd_ratio",
                         "n_above_99th_percentile")
} else {
  tier_level <- data.frame(
    tier = c("1. Autophagosome donor","2. Granule donor",
             "3. Shared lysosomal acceptor","4. Upstream signalling"),
    n_posthoc_genes = c(51, 16, 25, 15),
    median_sd_ratio = c(1.19, 1.04, 1.03, 0.91),
    n_above_99th_percentile = c(1, 0, 0, 0))
}

addWorksheet(wb, "S3_summaries")
writeData(wb, "S3_summaries",
  paste("Supplementary Table S3. Set-level and tier-level summaries.",
        "Panel A: pre-registered sets against 10,000 size-matched random gene sets",
        "drawn from the genome-wide background. Panel B: exploratory tiers,",
        "post-hoc genes only; genome-wide median 1.077, 99th percentile 1.882."),
  startRow = 1)
addStyle(wb, "S3_summaries",
         createStyle(fontSize = 9, textDecoration = "italic", wrapText = TRUE),
         rows = 1, cols = 1)
setRowHeights(wb, "S3_summaries", rows = 1, heights = 40)

writeData(wb, "S3_summaries", "A. Pre-registered set-level tests",
          startRow = 3)
addStyle(wb, "S3_summaries", createStyle(textDecoration = "bold"),
         rows = 3, cols = 1)
writeData(wb, "S3_summaries", rnd(set_level), startRow = 4,
          headerStyle = hdr)

r0 <- 4 + nrow(set_level) + 3
writeData(wb, "S3_summaries", "B. Exploratory tiers (post-hoc, null result)",
          startRow = r0 - 1)
addStyle(wb, "S3_summaries", createStyle(textDecoration = "bold"),
         rows = r0 - 1, cols = 1)
writeData(wb, "S3_summaries", rnd(tier_level), startRow = r0,
          headerStyle = hdr)
setColWidths(wb, "S3_summaries", cols = 1:8, widths = "auto")


# ---- S4: post-hoc comparator panels ----------------------------------

bits <- list()

f <- file.path(FIGS, "exocytosis_panel.csv")
if (file.exists(f)) {
  x <- read.csv(f, stringsAsFactors = FALSE)
  x$panel <- "Regulated exocytosis"
  bits[[length(bits) + 1]] <- x[, intersect(
    c("panel","gene","group","status","sd_ratio","ci_lo","ci_hi","p","pct"),
    names(x))]
}

f <- file.path(FIGS, "rab_regulatory_layer.csv")
if (file.exists(f)) {
  x <- read.csv(f, stringsAsFactors = FALSE)
  names(x)[names(x) == "label"] <- "gene"
  names(x)[names(x) == "layer"] <- "group"
  x$panel <- "Rab GEF / GAP layer"
  bits[[length(bits) + 1]] <- x[, intersect(
    c("panel","gene","group","status","sd_ratio","ci_lo","ci_hi","p","pct"),
    names(x))]
}

f <- file.path(FIGS, "munc13_paralogues.csv")
if (file.exists(f)) {
  x <- read.csv(f, stringsAsFactors = FALSE)
  names(x)[names(x) == "label"] <- "gene"
  x$panel <- "Munc13 paralogues"
  x$group <- "Munc13 family"
  bits[[length(bits) + 1]] <- x[, intersect(
    c("panel","gene","group","sd_ratio","ci_lo","ci_hi","p","pct"),
    names(x))]
}

if (length(bits)) {
  s4 <- do.call(rbind, lapply(bits, function(d) {
    miss <- setdiff(c("panel","gene","group","status","sd_ratio",
                      "ci_lo","ci_hi","p","pct"), names(d))
    for (m in miss) d[[m]] <- NA
    d[, c("panel","gene","group","status","sd_ratio","ci_lo","ci_hi",
          "p","pct")]
  }))
  s4$gene <- gsub("\n", " ", s4$gene)

  add_sheet("S4_posthoc_panels",
    paste("Supplementary Table S4. Post-hoc comparator panels, all under the",
          "fully adjusted model. The regulated exocytosis panel tests whether",
          "elevated variance is a general property of beta-cell membrane fusion",
          "machinery (it is not: no gene above the 99th percentile). The Rab",
          "GEF/GAP panel is uninformative about Rab activity, since these",
          "proteins act catalytically and abundance is a poor proxy for",
          "nucleotide state. The Munc13 paralogue comparison is grouped by",
          "sequence homology rather than by any partition imposed after seeing",
          "results."),
    rnd(s4))
} else {
  message("No post-hoc panel CSVs found; S4 skipped.")
}


# ---- Write -----------------------------------------------------------

out <- file.path(OUT, "Supplementary_Tables.xlsx")
saveWorkbook(wb, out, overwrite = TRUE)

cat("Wrote", out, "\n")
cat("Sheets:", paste(names(wb), collapse = ", "), "\n")
cat("\nEach sheet carries its own explanatory note in row 1, so the\n")
cat("tables stand alone if separated from the manuscript.\n")
