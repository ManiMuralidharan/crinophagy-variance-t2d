# =====================================================================
# 33. Set-level tests under the 7-gene beta-cell score
#
# Completes section 5 of script 32, which was skipped because the frozen
# gene-set file is Entrez-keyed while the replication matrix is Ensembl-
# keyed. This maps between them and re-runs both pre-registered set-level
# tests under Model 2 (8-gene score) and Model 2b (7-gene score).
#
# Tests, both pre-registered:
#   (1) median shift  - set median SD ratio vs 10,000 size-matched
#                       random sets drawn from the same background
#   (2) tail enrichment - members above the model-specific 99th
#                       percentile vs a binomial expectation of 1%
#
# PREREQUISITES: mat, pd2, bg2 in session (or session_gse164416.RData
# on disk), and bg2b_beta7.RData written by script 32.
#
# RUN:  source("33_setlevel_ins_sensitivity.R")
# =====================================================================

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))

if (!all(sapply(c("mat", "pd2", "bg2"), exists)))
  load("session_gse164416.RData")
if (!exists("bg2b")) load("bg2b_beta7.RData")
stopifnot(exists("bg2"), exists("bg2b"))

set.seed(1)
NPERM <- 10000


# ---- 1. Read the frozen sets ----------------------------------------

SETFILE <- "frozen_gene_sets_v2_clean.txt"
stopifnot(file.exists(SETFILE))

gs <- read.delim(SETFILE, stringsAsFactors = FALSE)
cat("Frozen set file columns:", paste(names(gs), collapse = ", "), "\n")
print(table(gs$set))

ent_col <- grep("ENTREZ", names(gs), ignore.case = TRUE, value = TRUE)[1]
sym_col <- grep("SYMBOL", names(gs), ignore.case = TRUE, value = TRUE)[1]
stopifnot(!is.na(ent_col))
cat("Using key column:", ent_col, "\n")


# ---- 2. Entrez -> Ensembl -------------------------------------------
# Mapped programmatically. The mapping is 1:many, so a gene can return
# several Ensembl accessions; we keep those actually present in the
# matrix, and where more than one survives we keep the one with the
# highest mean expression and record the ambiguity.

suppressPackageStartupMessages({
  library(org.Hs.eg.db); library(AnnotationDbi)
})
keys <- unique(as.character(gs[[ent_col]]))
m <- AnnotationDbi::select(org.Hs.eg.db, keys = keys,
                           keytype = "ENTREZID", columns = "ENSEMBL")
detach("package:org.Hs.eg.db", unload = TRUE)
detach("package:AnnotationDbi", unload = TRUE)

m <- m[!is.na(m$ENSEMBL) & m$ENSEMBL %in% rownames(mat), ]

mean_expr <- rowMeans(mat[unique(m$ENSEMBL), , drop = FALSE])
m$expr <- mean_expr[m$ENSEMBL]
m <- m[order(m$ENTREZID, -m$expr), ]

dups <- unique(m$ENTREZID[duplicated(m$ENTREZID)])
if (length(dups)) {
  cat("\nEntrez IDs with >1 Ensembl accession in the matrix (",
      length(dups), "). Highest-expressed retained:\n", sep = "")
  print(m[m$ENTREZID %in% dups, c("ENTREZID", "ENSEMBL", "expr")],
        row.names = FALSE, digits = 4)
}
m <- m[!duplicated(m$ENTREZID), ]

gs$ENSEMBL <- m$ENSEMBL[match(as.character(gs[[ent_col]]), m$ENTREZID)]

cat("\n===== MAPPING SUMMARY =====\n")
for (s in unique(gs$set)) {
  sub <- gs[gs$set == s, ]
  cat(sprintf("  %-16s %2d members, %2d mapped and present\n",
              s, nrow(sub), sum(!is.na(sub$ENSEMBL))))
  lost <- sub[is.na(sub$ENSEMBL), ]
  if (nrow(lost)) {
    lab <- if (!is.na(sym_col)) lost[[sym_col]] else lost[[ent_col]]
    cat("      absent:", paste(lab, collapse = ", "), "\n")
  }
}

sets <- lapply(split(gs$ENSEMBL, gs$set), function(v) intersect(v, names(bg2b)))


# ---- 3. The two set-level tests -------------------------------------

median_shift <- function(bg, members, nperm = NPERM) {
  obs  <- median(bg[members])
  null <- replicate(nperm, median(sample(bg, length(members))))
  data.frame(n = length(members),
             observed = obs,
             null_median = median(null),
             P = (sum(null >= obs) + 1) / (nperm + 1))
}

tail_enrich <- function(bg, members) {
  thr <- quantile(bg, .99)
  k   <- sum(bg[members] > thr)
  data.frame(n = length(members),
             threshold = unname(thr),
             observed = k,
             expected = 0.01 * length(members),
             P = binom.test(k, length(members), 0.01,
                            alternative = "greater")$p.value,
             genes = paste(names(sort(bg[members][bg[members] > thr],
                                      decreasing = TRUE)), collapse = ", "))
}

run_all <- function(bg, label) {
  cat("\n=====", label, "=====\n")
  cat("-- median shift --\n")
  ms <- do.call(rbind, lapply(names(sets), function(s)
    cbind(set = s, median_shift(bg, sets[[s]]))))
  print(ms, row.names = FALSE, digits = 4)

  cat("-- tail enrichment (99th percentile) --\n")
  te <- do.call(rbind, lapply(names(sets), function(s)
    cbind(set = s, tail_enrich(bg, sets[[s]]))))
  print(te[, setdiff(names(te), "genes")], row.names = FALSE, digits = 4)
  for (i in seq_len(nrow(te)))
    cat(sprintf("     %s above threshold: %s\n", te$set[i],
                ifelse(nzchar(te$genes[i]), te$genes[i], "none")))

  list(median = ms, tail = te, model = label)
}

r8 <- run_all(bg2,  "MODEL 2  (8-gene beta score, primary)")
r7 <- run_all(bg2b, "MODEL 2b (7-gene beta score, INS removed)")


# ---- 4. Side by side ------------------------------------------------

cmp <- merge(
  data.frame(set = r8$median$set, P_median_8 = r8$median$P,
             med_8 = r8$median$observed,
             tail_8 = r8$tail$observed, P_tail_8 = r8$tail$P),
  data.frame(set = r7$median$set, P_median_7 = r7$median$P,
             med_7 = r7$median$observed,
             tail_7 = r7$tail$observed, P_tail_7 = r7$tail$P),
  by = "set")

cat("\n===== SIDE BY SIDE =====\n")
print(cmp, row.names = FALSE, digits = 4)

write.csv(cmp, "setlevel_ins_sensitivity.csv", row.names = FALSE)
write.csv(gs[, c(ent_col, sym_col, "set", "ENSEMBL")],
          "frozen_sets_entrez_to_ensembl.csv", row.names = FALSE)

cat("\nWrote setlevel_ins_sensitivity.csv and",
    "frozen_sets_entrez_to_ensembl.csv\n")
cat("\nNOTE. Permutation P values move with the seed. This script uses\n")
cat("set.seed(1) and 10,000 permutations. Whatever appears in the\n")
cat("manuscript must come from a recorded seed and this many shuffles.\n")
cat("The crinophagy median-shift P under Model 2 printed above is the\n")
cat("number that belongs in section 2.4 and Figure 2C.\n")
