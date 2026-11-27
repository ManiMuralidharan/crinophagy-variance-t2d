# =====================================================================
# 35. Restore a clean session and finish the INS sensitivity numbers
#
# Run this in a FRESH R session. It loads everything from disk, rebuilds
# the gene sets from the v3 file, and prints the two numbers still
# needed for the Methods 4.6 paragraph:
#   - crinophagy median shift under Model 2b at 1e6 permutations
#   - the matching null median
#
# RUN:  source("35_restore_and_run.R")
# Takes a few minutes at the permutation step.
# =====================================================================

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))

cat("Files present:\n")
for (f in c("session_gse164416.RData", "bg2b_beta7.RData",
            "frozen_gene_sets_v3.txt"))
  cat(sprintf("  %-28s %s\n", f, ifelse(file.exists(f), "yes", "MISSING")))

stopifnot(file.exists("session_gse164416.RData"),
          file.exists("frozen_gene_sets_v3.txt"))


# ---- 1. Restore ------------------------------------------------------

load("session_gse164416.RData")          # mat, pd2, bg1, bg2, ratio_all

if (file.exists("bg2b_beta7.RData")) {
  load("bg2b_beta7.RData")               # bg2b
} else {
  cat("\nbg2b_beta7.RData missing. Rebuilding the 7-gene background.\n")
  suppressPackageStartupMessages({
    library(org.Hs.eg.db); library(AnnotationDbi) })
  b7 <- AnnotationDbi::select(org.Hs.eg.db,
          keys = c("IAPP","NKX6-1","PDX1","MAFA","G6PC2","ADCYAP1","SLC30A8"),
          keytype = "SYMBOL", columns = "ENSEMBL")
  detach("package:org.Hs.eg.db", unload = TRUE)
  detach("package:AnnotationDbi", unload = TRUE)
  b7 <- unique(b7$ENSEMBL[!is.na(b7$ENSEMBL) & b7$ENSEMBL %in% rownames(mat)])
  pd2$beta7 <- colMeans(mat[b7, rownames(pd2), drop = FALSE])
  X   <- model.matrix(~ batch + imm + beta7, data = pd2)
  res <- lm.fit(X, t(mat[, rownames(pd2), drop = FALSE]))$residuals
  bg2b <- apply(res[pd2$group == "T2D", , drop = FALSE], 2, sd) /
          apply(res[pd2$group == "ND",  , drop = FALSE], 2, sd)
  save(bg2b, file = "bg2b_beta7.RData")
  cat("Rebuilt and saved bg2b_beta7.RData\n")
}

cat("\nRestored:", nrow(pd2), "donors,", nrow(mat), "genes\n")
print(table(pd2$group))


# ---- 2. Rebuild the gene sets from v3 --------------------------------

gs3 <- read.delim("frozen_gene_sets_v3.txt", stringsAsFactors = FALSE)
cat("\nv3 file:\n"); print(table(gs3$set))

suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi) })
mm <- AnnotationDbi::select(org.Hs.eg.db,
        keys = unique(as.character(gs3$ENTREZID)),
        keytype = "ENTREZID", columns = "ENSEMBL")
detach("package:org.Hs.eg.db", unload = TRUE)
detach("package:AnnotationDbi", unload = TRUE)

mm <- mm[!is.na(mm$ENSEMBL) & mm$ENSEMBL %in% rownames(mat), ]
mm$expr <- rowMeans(mat[mm$ENSEMBL, , drop = FALSE])
mm <- mm[order(mm$ENTREZID, -mm$expr), ]
mm <- mm[!duplicated(mm$ENTREZID), ]
gs3$ENSEMBL <- mm$ENSEMBL[match(as.character(gs3$ENTREZID), mm$ENTREZID)]

sets  <- lapply(split(gs3$ENSEMBL, gs3$set), function(v) intersect(v, names(bg2b)))
crino <- sets[["crinophagy"]]
macro <- sets[["macroautophagy"]]

cat("\nSet sizes in matrix:\n"); print(sapply(sets, length))
if (length(crino) != 31 || length(macro) != 14)
  stop("Expected 31 and 14. STOP and report this before using any result.")


# ---- 3. The two numbers ----------------------------------------------

cat("\nRunning 1e6 permutations. This takes a few minutes.\n")
obs <- median(bg2b[crino])
set.seed(1)
null <- replicate(1e6, median(sample(bg2b, length(crino))))

cat("\n===== MODEL 2b (7-gene beta score), crinophagy set =====\n")
cat(sprintf("  observed median SD ratio : %.3f\n", obs))
cat(sprintf("  null median              : %.3f\n", median(null)))
cat(sprintf("  permutation P            : %.4f\n",
            (sum(null >= obs) + 1) / (1e6 + 1)))

thr <- quantile(bg2b, .99)
above <- names(sort(bg2b[crino][bg2b[crino] > thr], decreasing = TRUE))
cat(sprintf("\n  99th percentile threshold: %.3f\n", thr))
cat(sprintf("  members above threshold  : %d of %d\n",
            length(above), length(crino)))
cat("  binomial P               :",
    signif(binom.test(length(above), length(crino), 0.01,
                      alternative = "greater")$p.value, 3), "\n")

save(mat, pd2, bg1, bg2, bg2b, sets, crino, macro,
     file = "session_full.RData")
cat("\nSession saved as session_full.RData.\n")
cat("Next time: load('session_full.RData') restores everything at once.\n")
