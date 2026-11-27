# =====================================================================
# 60_fusion_noncircular.R
# The fusion grouping is circular: UNC13D and STX7 are in it because
# they ranked 1st and 2nd. This removes them and asks whether the REST
# of the fusion machinery is still elevated. That is the only version
# of the test that is not question-begging.
# RUN:  source("60_fusion_noncircular.R")
# =====================================================================
setwd(Sys.getenv("HUMAN_ISLETS_DIR", ".")); set.seed(1)
load("session_full.RData")

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
ens <- function(s) { g <- unname(s2e[s]); g[!is.na(g) & g %in% names(bg2)] }

tst <- function(sym, lab, B = 1e5) {
  g <- ens(sym); obs <- median(bg2[g])
  nul <- replicate(B, median(sample(bg2, length(g))))
  cat(sprintf("  %-44s n=%2d  median %.4f  P = %.5f\n",
              lab, length(g), obs, (sum(nul >= obs) + 1)/(B + 1)))
  invisible(obs)
}

TOP2 <- c("UNC13D","STX7")
cat("\n===== THE NON-CIRCULAR TEST =====\n")
tst(setdiff(FUSION, TOP2),            "Fusion WITHOUT UNC13D and STX7 (8)")
tst(setdiff(cr$SYMBOL, FUSION),       "Non-fusion remainder (21)")
tst(setdiff(cr$SYMBOL, TOP2),         "Whole set without UNC13D and STX7 (29)")
cat("\n  If the first row is still clearly elevated and well below the third,\n")
cat("  the concentration survives removing the two genes that defined the\n")
cat("  grouping. If it sits near the third, it does not.\n")

cat("\n===== every fusion member, so the reader can see the spread =====\n")
fm <- data.frame(gene = FUSION,
                 sd_ratio = sapply(FUSION, function(s) {
                   g <- ens(s); if (!length(g)) NA else unname(bg2[g]) }),
                 pct = sapply(FUSION, function(s) {
                   g <- ens(s); if (!length(g)) NA else 100*mean(bg2 <= bg2[g]) }))
print(fm[order(-fm$sd_ratio), ], row.names = FALSE, digits = 4)
cat(sprintf("\n  members above the genome-wide median of %.3f: %d of %d\n",
            median(bg2), sum(fm$sd_ratio > median(bg2), na.rm = TRUE), nrow(fm)))

cat("\n===== stability: drop each member in turn =====\n")
loo <- sapply(FUSION, function(s) median(bg2[ens(setdiff(FUSION, s))]))
print(round(sort(loo), 4))
cat(sprintf("\n  full fusion median %.4f; range when one member is dropped %.4f to %.4f\n",
            median(bg2[ens(FUSION)]), min(loo), max(loo)))
cat("  A wide range means the median is carried by a few members, which at\n")
cat("  n = 10 it usually is. Report this rather than a single P value.\n")
