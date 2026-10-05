# =====================================================================
# Step 2: Modelling the crinophagy signature
#
# Cohort: 129 HumanIslets donors with proteomics (all T1_diabetes = 0)
#   diagnosis_computed: None 91, Pre.T2D 20, Type2 18
#   HbA1c available in 117
#
# Primary test: continuous HbA1c (better powered than 3-group)
# Secondary:    variance heterogeneity by group
# Control:      identical pipeline on macroautophagy set
# Anchor:       score vs insulin content per IEQ
# =====================================================================

library(GSVA)
library(dplyr)
library(ggplot2)
library(car)        # leveneTest
library(broom)

setwd("/Users/mmani3/Downloads/Human Islets")

prot_proc  <- readRDS("./proteomics/proc_prot.rds")
donor_meta <- read.csv("donor.csv")
gene_sets  <- read.delim("frozen_gene_sets_v2_clean.txt")


# ---- 1. Build the expression matrix ---------------------------------

expr <- as.matrix(prot_proc[, -1, drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- as.character(prot_proc[[1]])   # Entrez IDs

# proc_prot is Entrez-keyed, so use ENTREZID from the frozen file
sets <- split(as.character(gene_sets$ENTREZID), gene_sets$set)

cat("Set sizes in matrix:\n")
print(sapply(sets, function(s) sum(s %in% rownames(expr))))


# ---- 2. ssGSEA ------------------------------------------------------
# NOTE: ssgsea.norm = FALSE. The earlier run compressed all scores into
# 0.84-0.97, which is what the default normalisation does to a small
# gene set. Unnormalised scores retain dynamic range. Both are reported.

par_norm <- ssgseaParam(expr, sets, normalize = TRUE)
par_raw  <- ssgseaParam(expr, sets, normalize = FALSE)

sc_norm <- t(gsva(par_norm))
sc_raw  <- t(gsva(par_raw))

scores <- data.frame(
  record_id      = rownames(sc_raw),
  crino_raw      = sc_raw[, "crinophagy"],
  macro_raw      = sc_raw[, "macroautophagy"],
  crino_norm     = sc_norm[, "crinophagy"],
  macro_norm     = sc_norm[, "macroautophagy"],
  stringsAsFactors = FALSE
)

cat("\nScore ranges (raw):\n")
print(summary(scores[, c("crino_raw", "macro_raw")]))


# ---- 3. Merge with metadata -----------------------------------------

d <- scores %>%
  inner_join(donor_meta, by = "record_id") %>%
  mutate(
    diagnosis_computed = factor(diagnosis_computed,
                                levels = c("None", "Pre.T2D", "Type2")),
    donorsex = factor(donorsex)
  )

stopifnot(nrow(d) == 129)


# ---- 4. PRIMARY: continuous HbA1c -----------------------------------
# Better powered than the 3-group test (117 donors, continuous exposure)

covars <- "donorage + donorsex + bodymassindex + predistributionculturetime + percentieqrecoverypostculture"

m_crino <- lm(as.formula(paste("crino_raw ~ hba1c +", covars)), data = d)
m_macro <- lm(as.formula(paste("macro_raw ~ hba1c +", covars)), data = d)

cat("\n===== PRIMARY: crinophagy ~ HbA1c =====\n")
print(tidy(m_crino) %>% filter(term == "hba1c"))
cat("n =", nobs(m_crino), "\n")

cat("\n===== CONTROL: macroautophagy ~ HbA1c =====\n")
print(tidy(m_macro) %>% filter(term == "hba1c"))

# The specificity claim lives here: crinophagy term significant,
# macroautophagy term not. If both move, the effect is generic.


# ---- 5. SECONDARY: group means (the null result, reported honestly) --

m_grp <- lm(as.formula(paste("crino_raw ~ diagnosis_computed +", covars)),
            data = d)
cat("\n===== Group comparison (adjusted) =====\n")
print(anova(m_grp))


# ---- 6. VARIANCE HETEROGENEITY --------------------------------------
# The boxplot suggested T2D is more dispersed. This tests it formally,
# on residuals so that covariate effects are removed first.

d_resid <- d %>%
  filter(!is.na(crino_raw), !is.na(donorage), !is.na(bodymassindex),
         !is.na(predistributionculturetime),
         !is.na(percentieqrecoverypostculture))

m_adj <- lm(as.formula(paste("crino_raw ~", covars)), data = d_resid)
d_resid$resid_crino <- residuals(m_adj)

cat("\n===== VARIANCE BY GROUP (adjusted residuals) =====\n")
print(d_resid %>%
        group_by(diagnosis_computed) %>%
        summarise(n = n(),
                  sd = sd(resid_crino),
                  iqr = IQR(resid_crino)))

cat("\nLevene (median-centred, i.e. Brown-Forsythe):\n")
print(leveneTest(resid_crino ~ diagnosis_computed, data = d_resid))

# Same test on the comparator, for specificity
m_adj_mac <- lm(as.formula(paste("macro_raw ~", covars)), data = d_resid)
d_resid$resid_macro <- residuals(m_adj_mac)
cat("\nComparator variance test:\n")
print(leveneTest(resid_macro ~ diagnosis_computed, data = d_resid))


# ---- 7. OUTLIER: the 0.56 donor -------------------------------------

out <- d %>%
  mutate(z = as.numeric(scale(crino_raw))) %>%
  filter(abs(z) > 3) %>%
  select(record_id, crino_raw, z, diagnosis_computed, hba1c,
         predistributionculturetime, percentieqrecoverypostculture,
         pdinsulinperieq)

cat("\n===== OUTLIERS (|z| > 3) =====\n")
print(out)
cat("Check culture time and IEQ recovery before excluding.\n")
cat("Rerun sections 4 and 6 without them as a sensitivity analysis.\n")


# ---- 8. FUNCTIONAL ANCHORING ----------------------------------------
# This is the strongest available claim: does the signature track a
# physiological readout in the same donors?

for (fn in c("pdinsulinperieq", "pdinsulindnaratio", "pdinsulincontent")) {
  f <- as.formula(paste("crino_raw ~", fn, "+", covars))
  m <- lm(f, data = d)
  cat("\n----- crinophagy ~", fn, "-----\n")
  print(tidy(m) %>% filter(term == fn))
}


# ---- 9. SENSITIVITY: drop UNC13D ------------------------------------
# UNC13D was detected in 75.2% of donors. Confirm nothing hinges on it.

sets_no_unc <- sets
sets_no_unc$crinophagy <- setdiff(sets$crinophagy, "201294")

sc_no_unc <- t(gsva(ssgseaParam(expr, sets_no_unc, normalize = FALSE)))
d$crino_no_unc <- sc_no_unc[d$record_id, "crinophagy"]

m_sens <- lm(as.formula(paste("crino_no_unc ~ hba1c +", covars)), data = d)
cat("\n===== SENSITIVITY: UNC13D removed =====\n")
print(tidy(m_sens) %>% filter(term == "hba1c"))
cat("Correlation with full score:",
    round(cor(d$crino_raw, d$crino_no_unc, use = "complete.obs"), 3), "\n")


# ---- 10. Figures ----------------------------------------------------

p1 <- ggplot(d, aes(x = hba1c, y = crino_raw)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", colour = "#2e7d32") +
  labs(x = "HbA1c (%)", y = "Crinophagy ssGSEA score",
       title = "Crinophagy signature across the glycemic continuum") +
  theme_minimal(base_size = 12)

p2 <- d_resid %>%
  ggplot(aes(x = diagnosis_computed, y = resid_crino,
             fill = diagnosis_computed)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.6) +
  labs(x = NULL, y = "Adjusted crinophagy score (residual)",
       title = "Dispersion, not location, differs by glycemic status") +
  theme_minimal(base_size = 12) + theme(legend.position = "none")

ggsave("fig1_hba1c_continuous.pdf", p1, width = 6, height = 5)
ggsave("fig2_variance_by_group.pdf", p2, width = 6, height = 5)

write.csv(d, "modelling_dataset.csv", row.names = FALSE)
