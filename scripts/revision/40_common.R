# =====================================================================
# 40_common.R   Shared setup and helpers for the revision analyses
#
# Run this first in every session:   source("40_common.R")
# It loads nothing you have not already produced and defines the
# helpers used by scripts 41 to 44.
# =====================================================================

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))
set.seed(1)

suppressPackageStartupMessages({
  library(car)
})

## ---- 1. Restore -----------------------------------------------------

stopifnot(file.exists("session_full.RData"))
load("session_full.RData")          # mat, pd2, bg1, bg2, bg2b, sets, crino, macro

need <- c("mat", "pd2", "crino", "macro")
miss <- need[!sapply(need, exists)]
if (length(miss))
  stop("Missing objects: ", paste(miss, collapse = ", "),
       "\nRe-run 35_restore_and_run.R first.")

SAMPLES <- rownames(pd2)
stopifnot(all(SAMPLES %in% colnames(mat)))
stopifnot(all(c("group", "batch", "imm", "beta") %in% names(pd2)))

cat(sprintf("Restored: %d donors, %d genes\n", nrow(pd2), nrow(mat)))
print(table(pd2$group, pd2$batch))


## ---- 2. Residuals under a given model -------------------------------
# Fitted ONCE. Because no model contains `group`, permuting the group
# labels does not change these residuals. That is what makes the
# label-permutation null below cheap.

residual_matrix <- function(formula) {
  X <- model.matrix(formula, data = pd2)
  Y <- t(mat[, SAMPLES, drop = FALSE])          # donors x genes
  stopifnot(nrow(X) == nrow(Y))
  lm.fit(X, Y)$residuals                        # donors x genes
}

M1 <- ~ batch                                   # batch only
M2 <- ~ batch + imm + beta                      # fully adjusted (primary)

RES1 <- residual_matrix(M1)
RES2 <- residual_matrix(M2)


## ---- 3. SD ratio for every gene, given a labelling ------------------

sd_ratios <- function(RES, lab) {
  a <- RES[lab == "T2D", , drop = FALSE]
  b <- RES[lab == "ND",  , drop = FALSE]
  matrixStats_sd <- function(M) sqrt(colSums(sweep(M, 2, colMeans(M))^2) / (nrow(M) - 1))
  matrixStats_sd(a) / matrixStats_sd(b)
}

OBS_LAB <- as.character(pd2$group)
stopifnot(all(sort(unique(OBS_LAB)) == c("ND", "T2D")))


## ---- 4. Within-batch label permutation ------------------------------
# Library preparation batch is strongly associated with glycemic status
# (chi-square P = 0.00068), so labels are exchangeable only WITHIN batch.
# Unrestricted permutation would be the wrong null and a referee would
# say so.

BATCH_IDX <- split(seq_along(OBS_LAB), pd2$batch)

permute_labels <- function() {
  out <- OBS_LAB
  for (ix in BATCH_IDX) if (length(ix) > 1) out[ix] <- sample(OBS_LAB[ix])
  out
}

# Batches containing only one group contribute no permutation entropy.
# Report how many donors are actually permutable, because that is the
# effective basis of the test and the referee will ask.
informative <- sapply(BATCH_IDX, function(ix) length(unique(OBS_LAB[ix])) > 1)
cat(sprintf("\nBatches with both groups: %d of %d; permutable donors: %d of %d\n",
            sum(informative), length(BATCH_IDX),
            sum(lengths(BATCH_IDX[informative])), length(OBS_LAB)))
if (sum(informative) == 0)
  stop("No batch contains both groups. Within-batch permutation is not possible; ",
       "report this and use an alternative design.")

cat("\n40_common.R loaded. RES1, RES2, sd_ratios(), permute_labels() available.\n")
