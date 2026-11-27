# =====================================================================
# 43_multiplicity.R   C-9 and C-16
#
# C-16 The exploratory FDR family currently contains the 14
#      prospectively defined comparator genes. Methods 4.10 says frozen
#      members were excluded; 29 crinophagy genes were, the comparator
#      genes were not. Recompute across the 93 genuinely post-hoc genes.
#
# C-9  Corrected and uncorrected evidence for every panel gene under
#      both models, in one table, so a reader can apply either standard.
#
# This script reads your existing Supplementary_Tables.xlsx. It changes
# nothing in that file; it writes new CSVs beside it.
#
# RUN:  source("43_multiplicity.R")      (does not need 40_common.R)
# =====================================================================

setwd(Sys.getenv("HUMAN_ISLETS_DIR", "."))

if (!requireNamespace("readxl", quietly = TRUE))
  stop("install.packages('readxl') first")
library(readxl)

XL <- "./ijms-4686697-supplementary/Supplementary_Tables.xlsx"
stopifnot(file.exists(XL))

S1 <- as.data.frame(read_excel(XL, sheet = "S1_pre-registered_sets", skip = 2))
S2 <- as.data.frame(read_excel(XL, sheet = "S2_exploratory_panel",  skip = 2))
cat(sprintf("S1: %d rows   S2: %d rows\n", nrow(S1), nrow(S2)))


## ---- C-16: rebuild the exploratory family ---------------------------

comparator <- S1$gene[S1$set == "comparator"]
crinophagy <- S1$gene[S1$set == "crinophagy"]
cat(sprintf("\nFrom S1: %d crinophagy, %d comparator genes\n",
            length(crinophagy), length(comparator)))

cat("\nBEFORE (as submitted):\n")
print(table(S2$frozen_set_status))

in_comp <- S2$gene %in% comparator
in_cri  <- S2$gene %in% crinophagy
cat(sprintf("\n  of the %d S2 genes labelled post-hoc, %d are prospectively\n",
            sum(S2$frozen_set_status != "Frozen crinophagy set"),
            sum(in_comp & S2$frozen_set_status != "Frozen crinophagy set")))
cat("  defined comparator members. THIS IS THE ERROR.\n")

# corrected provenance label
S2$provenance <- ifelse(in_cri,  "Frozen crinophagy set",
                 ifelse(in_comp, "Frozen macroautophagy comparator",
                                 "Post-hoc (added after the primary analysis)"))
cat("\nAFTER (corrected provenance):\n")
print(table(S2$provenance))

posthoc <- S2[S2$provenance == "Post-hoc (added after the primary analysis)", ]
cat(sprintf("\nGenuinely post-hoc genes: %d  (manuscript currently says 107)\n",
            nrow(posthoc)))

pcol <- if ("p" %in% names(S2)) "p" else stop("no 'p' column found in S2")
posthoc$q_corrected <- p.adjust(posthoc[[pcol]], method = "BH")

old_q <- suppressWarnings(min(as.numeric(S2$q), na.rm = TRUE))
new_q <- min(posthoc$q_corrected, na.rm = TRUE)
cat(sprintf("\n  minimum q as published (107-gene family): %.4f\n", old_q))
cat(sprintf("  minimum q recomputed    (%3d-gene family): %.4f\n",
            nrow(posthoc), new_q))
if (new_q < 0.05) {
  cat("  *** A gene now survives correction. STOP and read the result before\n",
      "  *** writing anything; this changes the Results, not only the Methods.\n")
} else {
  cat("  Conclusion unchanged: still no gene surviving correction.\n")
}

write.csv(S2[, c("gene", "subgroup", "tier", "provenance", pcol)],
          "C16_S2_corrected_provenance.csv", row.names = FALSE)
write.csv(posthoc, "C16_exploratory_FDR_93.csv", row.names = FALSE)


## ---- C-9: one multiplicity table for both panels --------------------

fam <- function(df, set_name, model) {
  d <- df[df$set == set_name, ]
  p <- as.numeric(d[[paste0("p_", model)]])
  data.frame(gene = d$gene, set = set_name, model = model,
             sd_ratio = as.numeric(d[[paste0("sd_ratio_", model)]]),
             ci_lo = as.numeric(d[[paste0("ci_lo_", model)]]),
             ci_hi = as.numeric(d[[paste0("ci_hi_", model)]]),
             pct   = as.numeric(d[[paste0("pct_", model)]]),
             p_nominal = p,
             q_within_own_panel = p.adjust(p, method = "BH"))
}

C9 <- do.call(rbind, lapply(c("M1", "M2"), function(m)
        rbind(fam(S1, "crinophagy", m), fam(S1, "comparator", m))))

# the joint family a sceptical reader would construct
for (m in c("M1", "M2")) {
  ix <- C9$model == m
  C9$q_joint_45[ix] <- p.adjust(C9$p_nominal[ix], method = "BH")
}

C9 <- C9[order(C9$model, C9$p_nominal), ]
write.csv(C9, "C9_multiplicity_both_panels.csv", row.names = FALSE)

cat("\n===== C-9  top of the joint ranking, Model 2 =====\n")
print(head(C9[C9$model == "M2", c("gene","set","sd_ratio","p_nominal",
                                  "q_within_own_panel","q_joint_45")], 10),
      row.names = FALSE, digits = 4)

cat("\n  Report q_within_own_panel AND q_joint_45. The GABARAP problem\n")
cat("  Referee 2 raised is visible in this table, which is the point:\n")
cat("  show it rather than let the referee find it.\n")
cat("\nWrote C9_multiplicity_both_panels.csv, C16_exploratory_FDR_93.csv,\n")
cat("C16_S2_corrected_provenance.csv\n")
