# Core-fusion sensitivity using the contemporaneous Fusion group from
# 15_volcano_and_df_check.R, not a post-hoc reconstruction.
setwd(Sys.getenv("HUMAN_ISLETS_DIR", ".")); set.seed(1)
load("session_full.RData")   # mat, pd2, bg2, crino, macro

FUSION <- c("UNC13D","STX7","STX8","STX17","VAMP7","VAMP8","SNAP23",
            "VPS33A","VPS16","VPS11")
gs <- read.delim("frozen_gene_sets_v3.txt", stringsAsFactors = FALSE)
cr <- gs[gs$set == "crinophagy", ]

suppressPackageStartupMessages({library(org.Hs.eg.db); library(AnnotationDbi)})
mm <- AnnotationDbi::select(org.Hs.eg.db, keys = as.character(cr$ENTREZID),
                            keytype = "ENTREZID", columns = "ENSEMBL")
detach("package:org.Hs.eg.db", unload = TRUE); detach("package:AnnotationDbi", unload = TRUE)
mm <- mm[!is.na(mm$ENSEMBL) & mm$ENSEMBL %in% names(bg2), ]
mm$expr <- rowMeans(mat[mm$ENSEMBL, , drop = FALSE])
mm <- mm[order(mm$ENTREZID, -mm$expr), ]; mm <- mm[!duplicated(mm$ENTREZID), ]
cr$ENSEMBL <- mm$ENSEMBL[match(as.character(cr$ENTREZID), mm$ENTREZID)]
s2e <- setNames(cr$ENSEMBL, cr$SYMBOL)

tst <- function(sym, lab, B = 1e5) {
  g <- intersect(unname(s2e[sym]), names(bg2)); g <- g[!is.na(g)]
  obs <- median(bg2[g]); nul <- replicate(B, median(sample(bg2, length(g))))
  cat(sprintf("  %-40s n=%2d  median %.4f  null %.4f  P = %.5f\n",
              lab, length(g), obs, median(nul), (sum(nul >= obs) + 1)/(B + 1)))
}
cat("\n== Model 2, random-gene-set null ==\n")
tst(cr$SYMBOL,                                   "Full frozen set (31)")
tst(FUSION,                                      "Fusion as defined (10)")
tst(setdiff(FUSION, "STX17"),                    "Fusion without STX17 (9)")
tst(setdiff(FUSION, "VAMP8"),                    "Fusion without VAMP8 (9)")
tst(setdiff(FUSION, c("STX17","VAMP8")),         "Fusion without both (8)")
tst(setdiff(cr$SYMBOL, FUSION),                  "Non-fusion remainder (21)")

cat("\nREAD THIS BEFORE QUOTING ANY OF IT:\n")
cat("  1. 13_canonical_analysis.R records that this partition was suggested\n")
cat("     by the results and that testing it is circular, because UNC13D and\n")
cat("     STX7 are in the fusion group BECAUSE they ranked first and second.\n")
cat("  2. The STX17 and VAMP8 variants answer the reviewer's question about\n")
cat("     set composition. If dropping them improves the P value, that is not\n")
cat("     a finding. Report all six rows together or none of them.\n")
cat("  3. These are random-gene-set P values: extremeness among comparable\n")
cat("     sets, not disease association.\n")
