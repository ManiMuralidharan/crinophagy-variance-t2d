# =====================================================================
# 52_cross_cohort.R
#
# Quantifies the claim the revised paper would rest on:
#   "UNC13D is elevated in both cohorts, and almost nothing else is."
#
# Genome-wide agreement between the two layers is weak (r = 0.06), so a
# single concordant gene could be chance. This script measures how
# unlikely the concordance actually is, rather than asserting it.
#
# Needs:  discovery_sd_ratios_all.rds   (from 51_discovery_gabarap.R)
#         session_full.RData            (bg1, bg2, mat, crino, macro)
#
# RUN:  source("52_cross_cohort.R")
# =====================================================================

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))
set.seed(1)

stopifnot(file.exists("discovery_sd_ratios_all.rds"),
          file.exists("session_full.RData"))

DISC <- readRDS("discovery_sd_ratios_all.rds")      # Entrez-keyed
load("session_full.RData")                          # bg1, bg2, mat, crino, macro
REP  <- bg2                                         # Model 2, primary
cat(sprintf("discovery: %d proteins   replication: %d genes\n",
            length(DISC), length(REP)))


## ---- 1. Entrez -> Ensembl, 1:many resolved by expression ------------

suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
map <- AnnotationDbi::select(org.Hs.eg.db, keys = names(DISC),
                             keytype = "ENTREZID", columns = "ENSEMBL")
detach("package:org.Hs.eg.db", unload = TRUE)
detach("package:AnnotationDbi", unload = TRUE)

map <- map[!is.na(map$ENSEMBL) & map$ENSEMBL %in% names(REP), ]
amb <- sum(duplicated(map$ENTREZID))
map$expr <- rowMeans(mat[map$ENSEMBL, , drop = FALSE])
map <- map[order(map$ENTREZID, -map$expr), ]
map <- map[!duplicated(map$ENTREZID), ]             # keep highest-expressed
cat(sprintf("mapped %d proteins to replication genes (%d ambiguous, highest-expressed kept)\n",
            nrow(map), amb))

d <- DISC[map$ENTREZID]
r <- REP[map$ENSEMBL]
names(d) <- names(r) <- map$ENTREZID
ok <- is.finite(d) & is.finite(r) & d > 0 & r > 0
d <- d[ok]; r <- r[ok]
N <- length(d)
cat(sprintf("proteins measured in BOTH cohorts: %d\n", N))
if (N < 1000) stop("intersection too small; check the mapping")


## ---- 2. Agreement ---------------------------------------------------

cat("\n===== agreement between layers =====\n")
cat(sprintf("  Pearson  on log SD ratio : %.3f\n", cor(log(d), log(r))))
cat(sprintf("  Spearman on SD ratio     : %.3f\n", cor(d, r, method = "spearman")))
cat("  (published genome-wide figure was r = 0.06)\n")


## ---- 3. Concordant tails --------------------------------------------
# Thresholds are computed WITHIN the intersected set. Using genome-wide
# thresholds here would bias the count, because the intersection is not
# a random subset of either cohort.

tail_test <- function(q, label) {
  td <- quantile(d, q); tr <- quantile(r, q)
  hi_d <- d >= td; hi_r <- r >= tr
  both <- hi_d & hi_r
  obs  <- sum(both)
  exp  <- N * (1 - q)^2
  ft   <- fisher.test(table(hi_d, hi_r), alternative = "greater")
  # permutation preserving both marginals exactly
  nperm <- 20000
  nullc <- replicate(nperm, sum(hi_d & sample(hi_r)))
  pperm <- (sum(nullc >= obs) + 1) / (nperm + 1)
  cat(sprintf("\n-- %s (discovery >= %.3f, replication >= %.3f) --\n", label, td, tr))
  cat(sprintf("   genes above in discovery   : %d\n", sum(hi_d)))
  cat(sprintf("   genes above in replication : %d\n", sum(hi_r)))
  cat(sprintf("   above in BOTH              : %d   (expected by chance %.1f)\n", obs, exp))
  cat(sprintf("   enrichment                 : %.2f-fold\n", obs / exp))
  cat(sprintf("   Fisher one-sided P         : %.4g\n", ft$p.value))
  cat(sprintf("   permutation P (marginals fixed): %.4g\n", pperm))
  list(q = q, thr_d = td, thr_r = tr, n_both = obs, expected = exp,
       p_fisher = ft$p.value, p_perm = pperm, genes = names(d)[both])
}

t99 <- tail_test(0.99, "top 1%")
t95 <- tail_test(0.95, "top 5%")


## ---- 4. Who are they? -----------------------------------------------

suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
sym <- AnnotationDbi::mapIds(org.Hs.eg.db, keys = t99$genes,
                             keytype = "ENTREZID", column = "SYMBOL")
detach("package:org.Hs.eg.db", unload = TRUE)
detach("package:AnnotationDbi", unload = TRUE)

conc <- data.frame(entrez = t99$genes, symbol = unname(sym),
                   discovery = unname(d[t99$genes]),
                   replication = unname(r[t99$genes]),
                   pct_disc = 100 * sapply(d[t99$genes], function(v) mean(d <= v)),
                   pct_rep  = 100 * sapply(r[t99$genes], function(v) mean(r <= v)))
conc <- conc[order(-conc$discovery), ]
cat("\n===== genes above the 99th percentile in BOTH cohorts =====\n")
print(conc, row.names = FALSE, digits = 4)

UNC13D <- "201294"
cat("\nUNC13D in this intersection: ")
if (UNC13D %in% names(d)) {
  cat(sprintf("discovery %.3f (%.2f pct), replication %.3f (%.2f pct)%s\n",
              d[UNC13D], 100*mean(d <= d[UNC13D]),
              r[UNC13D], 100*mean(r <= r[UNC13D]),
              ifelse(UNC13D %in% t99$genes, "  -- IN the concordant set", "")))
} else cat("NOT PRESENT (check the mapping)\n")


## ---- 5. The panel genes, both cohorts -------------------------------

PANEL <- c(UNC13D = "201294", RAB3A = "5864", GABARAP = "11337",
           CTSD = "1509", GABARAPL1 = "23710", GABARAPL2 = "11345")
pan <- do.call(rbind, lapply(names(PANEL), function(s) {
  e <- PANEL[[s]]
  if (!e %in% names(d)) return(data.frame(gene = s, discovery = NA, pct_disc = NA,
                                          replication = NA, pct_rep = NA))
  data.frame(gene = s, discovery = unname(d[e]), pct_disc = 100*mean(d <= d[e]),
             replication = unname(r[e]), pct_rep = 100*mean(r <= r[e]))
}))
cat("\n===== panel genes in both cohorts =====\n")
print(pan, row.names = FALSE, digits = 4)

write.csv(conc, "C_concordant_tail_genes.csv", row.names = FALSE)
write.csv(pan,  "C_panel_both_cohorts.csv", row.names = FALSE)
save(t99, t95, d, r, file = "C_cross_cohort.RData")

cat("\n===== HOW TO READ THIS =====\n")
cat("  If n_both is small (a handful) and the permutation P is small, then\n")
cat("  'UNC13D is one of very few genes elevated in both cohorts' is a real,\n")
cat("  quantified claim and the revised paper can be built on it.\n")
cat("  If n_both is large, concordance is common and UNC13D being in it means\n")
cat("  little. Write whichever the number supports.\n")
cat("  Note the discovery values are from a reconstruction accurate to about\n")
cat("  1-3%. That is fine for ranks and tails; do not quote them as the\n")
cat("  published point estimates.\n")
cat("\nWrote C_concordant_tail_genes.csv, C_panel_both_cohorts.csv, C_cross_cohort.RData\n")
