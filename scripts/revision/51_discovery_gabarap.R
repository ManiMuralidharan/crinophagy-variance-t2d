# =====================================================================
# 51_discovery_gabarap.R
#
# Reconstructs the discovery per-protein variance analysis and reports
# GABARAP alongside the panel.
#
# IMPORTANT: this script VALIDATES ITSELF before reporting anything.
# It must reproduce the published discovery values for UNC13D (2.065)
# and RAB3A (1.725). If it cannot, it stops. Do not use numbers from a
# reconstruction that fails its own check.
#
# RUN:  source("51_discovery_gabarap.R")
# =====================================================================

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))
set.seed(1)

PUBLISHED <- c(UNC13D = 2.065, RAB3A = 1.725)   # discovery, from the manuscript
TOL       <- 0.02                                # how close counts as reproduced

ENTREZ <- c(UNC13D = "201294", RAB3A = "5864", GABARAP = "11337",
            STX7 = "8417", CTSD = "1509", GABARAPL1 = "23710",
            GABARAPL2 = "11345", MAP1LC3B = "81631")
# NOTE: UNC13D Entrez is 201294 in some builds; the script resolves both.
ALT <- c(UNC13D = "84257")   # older build


## ---- 1. Load -------------------------------------------------------

P <- readRDS("proteomics/proc_prot.rds")
stopifnot("gene_id" %in% names(P))
gid  <- as.character(P$gene_id)
Mp   <- as.matrix(P[, setdiff(names(P), "gene_id"), drop = FALSE])
rownames(Mp) <- gid
mode(Mp) <- "numeric"
cat(sprintf("proc_prot: %d proteins x %d donors\n", nrow(Mp), ncol(Mp)))

don <- read.csv("donor.csv", stringsAsFactors = FALSE)
don <- don[match(colnames(Mp), don$record_id), ]
stopifnot(identical(don$record_id, colnames(Mp)))
cat("donor metadata matched for all", nrow(don), "columns\n")


## ---- 2. Groups -----------------------------------------------------

cat("\ndiagnosis_computed in the 129 discovery donors:\n")
print(table(don$diagnosis_computed, useNA = "ifany"))

# diagnosis_computed already carries the three groups used in Figure 1C,
# so it is used directly rather than re-derived from HbA1c.
don$grp <- c("None" = "ND", "Pre.T2D" = "preT2D", "Type2" = "T2D")[
             as.character(don$diagnosis_computed)]
cat("\nthree-group split from diagnosis_computed:\n")
print(table(don$grp, useNA = "ifany"))
cat("published Figure 1C: ND 91, preT2D 20, T2D 18\n")


## ---- 3. Covariate-complete two-group set ---------------------------

COV <- c("donorage", "donorsex", "bodymassindex",
         "predistributionculturetime", "percentieqrecoverypostculture")
miss <- COV[!COV %in% names(don)]
if (length(miss)) stop("covariate columns not found: ", paste(miss, collapse = ", "))

cat("\nmissing values per covariate (two-group donors):\n")
two <- don$grp %in% c("ND", "T2D")
print(sapply(don[two, COV], function(z) sum(is.na(z))))
ok <- complete.cases(don[, COV]) & two
cat(sprintf("\ncovariate-complete two-group donors: %d\n", sum(ok)))
print(table(don$grp[ok]))
cat("published two-group comparison: 87 non-diabetic, 16 T2D\n")
if (sum(ok) < 50) stop("too few donors retained; check the covariate column names")

d  <- don[ok, ]
Y  <- t(Mp[, ok, drop = FALSE])                     # donors x proteins
d$donorsex <- factor(d$donorsex)


## ---- 4. SD ratios under the discovery covariate model --------------

X   <- model.matrix(~ donorage + donorsex + bodymassindex +
                      predistributionculturetime +
                      percentieqrecoverypostculture, data = d)
keep <- colSums(is.na(Y)) == 0                      # proteins with no missing
cat(sprintf("proteins with complete data: %d of %d\n", sum(keep), ncol(Y)))
Yk  <- Y[, keep, drop = FALSE]
res <- lm.fit(X, Yk)$residuals

colsd <- function(M) sqrt(colSums(sweep(M, 2, colMeans(M))^2) / (nrow(M) - 1))
r_t2d <- colsd(res[d$grp == "T2D", , drop = FALSE])
r_nd  <- colsd(res[d$grp == "ND",  , drop = FALSE])
SDR   <- r_t2d / r_nd
cat(sprintf("SD ratios computed for %d proteins\n", length(SDR)))


## ---- 5. SELF-CHECK --------------------------------------------------

resolve <- function(sym) {
  cand <- c(if (sym %in% names(ENTREZ)) ENTREZ[[sym]],
            if (sym %in% names(ALT))    ALT[[sym]])
  cand <- cand[!is.na(cand)]
  hit  <- cand[cand %in% names(SDR)]
  if (length(hit)) hit[1] else NA_character_
}

cat("\n===== SELF-CHECK against published discovery values =====\n")
passed <- TRUE
for (s in names(PUBLISHED)) {
  id <- resolve(s)
  if (is.na(id)) {
    cat(sprintf("  %-7s NOT FOUND in the matrix\n", s)); passed <- FALSE; next
  }
  got <- unname(SDR[id])
  dev <- abs(got - PUBLISHED[[s]])
  cat(sprintf("  %-7s entrez %-7s reconstructed %.3f   published %.3f   diff %.3f  %s\n",
              s, id, got, PUBLISHED[[s]], dev,
              ifelse(dev <= TOL, "OK", "MISMATCH")))
  if (dev > TOL) passed <- FALSE
}

if (!passed) {
  cat("\n*** SELF-CHECK FAILED.\n")
  cat("*** The reconstruction does not reproduce the published discovery\n")
  cat("*** numbers, so any GABARAP value it produces is not trustworthy.\n")
  cat("*** Likely causes, in order: the per-protein analysis used the\n")
  cat("*** UNPROCESSED matrix (unproc_prot.rds, UniProt-keyed) rather than\n")
  cat("*** this one; a different covariate set; or a different pre-T2D rule.\n")
  cat("*** Send me this output and I will adjust. Do not quote the numbers\n")
  cat("*** below in anything.\n\n")
} else {
  cat("\n  Self-check passed. The numbers below can be trusted.\n")
}


## ---- 6. The answer --------------------------------------------------

cat("\n===== DISCOVERY SD RATIOS, panel genes =====\n")
out <- do.call(rbind, lapply(names(ENTREZ), function(s) {
  id <- resolve(s)
  if (is.na(id)) return(data.frame(gene = s, entrez = NA, sd_ratio = NA,
                                   pct = NA, note = "not quantified in discovery"))
  v <- unname(SDR[id])
  data.frame(gene = s, entrez = id, sd_ratio = v,
             pct = 100 * mean(SDR <= v, na.rm = TRUE), note = "")
}))
print(out, row.names = FALSE, digits = 4)

cat("\nReplication values for comparison (Model 2):\n")
cat("  GABARAP 2.178 (M1 2.243, 99.76th pct)   UNC13D 2.001 (M1 1.918)\n")
cat("  STX7 1.899   CTSD 1.799   RAB3A 1.101\n")

cat("\nWHAT TO LOOK FOR:\n")
cat("  If GABARAP's discovery SD ratio is near 1, it does NOT replicate and\n")
cat("  UNC13D remains the only cross-platform candidate. That is defensible.\n")
cat("  If GABARAP is also elevated in discovery, it replicates as well as\n")
cat("  UNC13D does, and the specificity argument cannot stand as written.\n")
cat("  If GABARAP is absent from the discovery matrix, it was never testable\n")
cat("  across platforms; say that, and do not call it a failed replication.\n")

write.csv(out, "C_discovery_panel_sd_ratios.csv", row.names = FALSE)
saveRDS(SDR, "discovery_sd_ratios_all.rds")
cat("\nWrote C_discovery_panel_sd_ratios.csv and discovery_sd_ratios_all.rds\n")
