# =====================================================================
# 34. Two loose ends
#
#  (a) The frozen file has 13 comparator genes; the manuscript, the
#      Methods list and Figure 2B all say 14. Identify the missing one.
#  (b) Settle the Model 1 (batch-only) set-level permutation P, which
#      the manuscript reports as 0.0066 and a later rebuild gave 0.0086.
#
# RUN:  source("34_reconcile.R")
# =====================================================================

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))

if (!all(sapply(c("mat", "pd2", "bg1", "bg2"), exists)))
  load("session_gse164416.RData")
stopifnot(exists("bg1"), exists("bg2"))
if (!exists("bg2b")) load("bg2b_beta7.RData")

set.seed(1)
NPERM <- 10000


# ---- (a) Which comparator gene is missing? --------------------------

# The 14 named in Methods 4.1
METHODS_14 <- c("ATG3", "ATG5", "ATG7", "BECN1", "ULK1", "WIPI2",
                "MAP1LC3A", "MAP1LC3B", "GABARAP", "GABARAPL1",
                "GABARAPL2", "SQSTM1", "NBR1", "CALCOCO2")

gs <- read.delim("frozen_gene_sets_v2_clean.txt", stringsAsFactors = FALSE)
frozen_comp <- gs$SYMBOL[gs$set == "macroautophagy"]

cat("===== COMPARATOR SET RECONCILIATION =====\n")
cat("Methods 4.1 names:", length(METHODS_14), "\n")
cat("Frozen file has: ", length(frozen_comp), "\n\n")
cat("Frozen file contents:\n  ",
    paste(sort(frozen_comp), collapse = ", "), "\n\n")

cat("In Methods but NOT in the frozen file:\n  ",
    paste(setdiff(METHODS_14, frozen_comp), collapse = ", "), "\n")
cat("In the frozen file but NOT in Methods:\n  ",
    paste(setdiff(frozen_comp, METHODS_14), collapse = ", "), "\n\n")

cat("ACTION. Decide which list is authoritative, then check it against\n")
cat("the Zenodo deposit (10.5281/zenodo.23060343). If the deposit holds\n")
cat("13, the manuscript must say 13 everywhere and Figure 2B and 2D must\n")
cat("be rebuilt. If the deposit holds 14, this local file is a stale copy\n")
cat("and should be replaced rather than corrected.\n")

# Also look for any other frozen-set file lying around
others <- list.files(".", pattern = "frozen.*\\.(txt|tsv|csv)$",
                     recursive = TRUE)
if (length(others) > 1) {
  cat("\nOther frozen-set files present. Check whether one holds 14:\n")
  for (f in others) {
    d <- try(read.delim(f, stringsAsFactors = FALSE), silent = TRUE)
    if (!inherits(d, "try-error") && "set" %in% names(d))
      cat(sprintf("  %-45s %s\n", f,
                  paste(names(table(d$set)), table(d$set),
                        sep = "=", collapse = "  ")))
  }
}


# ---- (b) Model 1 set-level tests ------------------------------------
# Same two tests, same seed, same number of shuffles, under the
# batch-only background.

suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi)
})
m <- AnnotationDbi::select(org.Hs.eg.db,
                           keys = unique(as.character(gs$ENTREZID)),
                           keytype = "ENTREZID", columns = "ENSEMBL")
detach("package:org.Hs.eg.db", unload = TRUE)
detach("package:AnnotationDbi", unload = TRUE)

m <- m[!is.na(m$ENSEMBL) & m$ENSEMBL %in% rownames(mat), ]
m$expr <- rowMeans(mat[m$ENSEMBL, , drop = FALSE])
m <- m[order(m$ENTREZID, -m$expr), ]
m <- m[!duplicated(m$ENTREZID), ]
gs$ENSEMBL <- m$ENSEMBL[match(as.character(gs$ENTREZID), m$ENTREZID)]

sets <- lapply(split(gs$ENSEMBL, gs$set),
               function(v) intersect(v, names(bg1)))

median_shift <- function(bg, members, nperm = NPERM) {
  obs  <- median(bg[members])
  null <- replicate(nperm, median(sample(bg, length(members))))
  data.frame(n = length(members), observed = obs,
             null_median = median(null),
             P = (sum(null >= obs) + 1) / (nperm + 1))
}
tail_enrich <- function(bg, members) {
  thr <- quantile(bg, .99); k <- sum(bg[members] > thr)
  data.frame(n = length(members), threshold = unname(thr), observed = k,
             expected = 0.01 * length(members),
             P = binom.test(k, length(members), 0.01,
                            alternative = "greater")$p.value)
}

cat("\n\n===== MODEL 1 (batch only) =====\n")
cat("-- median shift --\n")
print(do.call(rbind, lapply(names(sets), function(s)
  cbind(set = s, median_shift(bg1, sets[[s]])))),
  row.names = FALSE, digits = 4)
cat("-- tail enrichment --\n")
print(do.call(rbind, lapply(names(sets), function(s)
  cbind(set = s, tail_enrich(bg1, sets[[s]])))),
  row.names = FALSE, digits = 4)

cat("\nThe crinophagy median-shift P above is the number for Figure 4C\n")
cat("and section 2.5, replacing both 0.0066 and 0.0086. Record the seed.\n")


# ---- (c) Stability of that P across seeds ---------------------------
# A permutation P near 0.006 is estimated from ~60 exceedances in 10,000
# shuffles, so it carries real Monte Carlo error. Quantify it once so you
# can state the P to a defensible number of digits.

cat("\n===== SEED STABILITY, crinophagy under Model 1 =====\n")
crino <- sets[["crinophagy"]]
obs   <- median(bg1[crino])
ps <- sapply(1:10, function(s) {
  set.seed(s)
  null <- replicate(NPERM, median(sample(bg1, length(crino))))
  (sum(null >= obs) + 1) / (NPERM + 1)
})
print(round(ps, 5))
cat(sprintf("  range %.5f to %.5f, mean %.5f\n",
            min(ps), max(ps), mean(ps)))
cat("  Report to two significant figures, or raise NPERM to 1e5\n")
cat("  and report the seed. Do not report four digits of a quantity\n")
cat("  that moves in the third.\n")
set.seed(1)
