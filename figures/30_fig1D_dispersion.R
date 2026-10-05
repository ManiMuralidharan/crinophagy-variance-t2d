# =====================================================================
# Figure 1D: dispersion of the adjusted pathway score
#
# SELF-CONTAINED. Needs nothing in memory. Reads only:
#   modelling_dataset.csv   (written by script 02)
# in the working folder. If that file is missing, see the note at the
# bottom of this script.
#
# Produces the primary discovery result as its own panel:
#   crinophagy    SD 62.9 (T2D) vs 39.0 (non-diabetic)
#                 SD ratio 1.62, 95% CI 0.98 to 2.24
#                 Brown-Forsythe P = 0.018, permutation P = 0.013
#   comparator    no corresponding difference (P = 0.62)
#
# Run:
#   setwd("/Users/mmani3/Downloads/Human Islets")
#   source("30_fig1D_dispersion.R")
# =====================================================================

library(ggplot2)
library(dplyr)
library(car)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

set.seed(1)
NPERM <- 10000
NBOOT <- 10000

f <- "modelling_dataset.csv"
if (!file.exists(f))
  stop("modelling_dataset.csv not found in ", getwd(),
       "\nSee the note at the end of this script.")

# NOTE: read.csv turns the string "None" into NA by default, which would
# silently delete all 91 non-diabetic donors. Read it this way only.
d <- read.csv(f, na.strings = c("NA", ""), stringsAsFactors = FALSE)
d$diagnosis_computed[is.na(d$diagnosis_computed)] <- "None"

d$group <- factor(d$diagnosis_computed,
                  levels = c("None", "Pre.T2D", "Type2"),
                  labels = c("Non-diabetic", "Pre-T2D", "T2D"))

COV <- c("donorage", "donorsex", "bodymassindex",
         "predistributionculturetime", "percentieqrecoverypostculture")

for (v in c("crino_raw", "macro_raw", setdiff(COV, "donorsex")))
  d[[v]] <- as.numeric(d[[v]])

ok <- complete.cases(d[, c("crino_raw", "macro_raw", COV)])
d2 <- droplevels(d[ok, ])
d2$donorsex <- factor(d2$donorsex)

cat("Donors after covariate filtering:", nrow(d2), "\n")
print(table(d2$group))


# ---- residualise both pathway scores ---------------------------------

fm <- function(y) as.formula(paste(y, "~", paste(COV, collapse = " + ")))
d2$r_crino <- residuals(lm(fm("crino_raw"), data = d2))
d2$r_macro <- residuals(lm(fm("macro_raw"), data = d2))


# ---- statistics ------------------------------------------------------

stats_for <- function(col) {
  lev <- leveneTest(as.formula(paste(col, "~ group")), data = d2)
  a <- d2[[col]][d2$group == "T2D"]
  b <- d2[[col]][d2$group == "Non-diabetic"]

  # permutation on the Brown-Forsythe statistic
  vals <- d2[[col]]; labs <- d2$group
  bf <- function(v, g) {
    gs <- split(v, g); gs <- gs[sapply(gs, length) > 1]
    car::leveneTest(unlist(gs),
                    factor(rep(names(gs), sapply(gs, length))),
                    center = median)$`F value`[1]
  }
  obs <- bf(vals, labs)
  perm <- labs
  null <- numeric(NPERM)
  for (i in seq_len(NPERM)) { perm <- sample(labs); null[i] <- bf(vals, perm) }
  p_perm <- (sum(null >= obs) + 1) / (NPERM + 1)

  boot <- replicate(NBOOT,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))

  list(sd_nd = sd(b), sd_t2d = sd(a),
       ratio = sd(a) / sd(b),
       lo = unname(quantile(boot, .025)),
       hi = unname(quantile(boot, .975)),
       p_bf = lev$`Pr(>F)`[1], p_perm = p_perm)
}

cat("\nComputing (10,000 permutations, 10,000 bootstraps each)...\n")
s_crino <- stats_for("r_crino")
s_macro <- stats_for("r_macro")

cat("\n===== CRINOPHAGY =====\n")
cat(sprintf("  SD  non-diabetic %.1f   T2D %.1f\n", s_crino$sd_nd, s_crino$sd_t2d))
cat(sprintf("  SD ratio %.3f, 95%% CI %.3f to %.3f\n",
            s_crino$ratio, s_crino$lo, s_crino$hi))
cat(sprintf("  Brown-Forsythe P = %.4f   permutation P = %.4f\n",
            s_crino$p_bf, s_crino$p_perm))

cat("\n===== MACROAUTOPHAGY COMPARATOR =====\n")
cat(sprintf("  SD  non-diabetic %.1f   T2D %.1f\n", s_macro$sd_nd, s_macro$sd_t2d))
cat(sprintf("  SD ratio %.3f\n", s_macro$ratio))
cat(sprintf("  Brown-Forsythe P = %.4f\n", s_macro$p_bf))

cat("\nThese should match the manuscript: crinophagy SD 62.9 vs 39.0,\n")
cat("Brown-Forsythe P = 0.018, permutation P = 0.013; comparator P = 0.62.\n")
cat("If they differ, the manuscript numbers need updating, not the code.\n")


# ---- figure ----------------------------------------------------------

long <- bind_rows(
  d2 %>% transmute(group, resid = r_crino, pathway = "Crinophagy"),
  d2 %>% transmute(group, resid = r_macro, pathway = "Macroautophagy comparator")
) %>%
  group_by(pathway) %>%
  mutate(z = (resid - mean(resid)) / sd(resid)) %>%
  ungroup() %>%
  mutate(pathway = factor(pathway,
           levels = c("Crinophagy", "Macroautophagy comparator")))

# per-group SD of the standardised values, for the annotation
sd_tab <- long %>% group_by(pathway, group) %>%
  summarise(sd = sd(z), .groups = "drop")

ann <- data.frame(
  pathway = factor(c("Crinophagy", "Macroautophagy comparator"),
                   levels = levels(long$pathway)),
  txt = c(sprintf("SD ratio %.2f (95%% CI %.2f\u2013%.2f)\nBrown-Forsythe P = %.3f\npermutation P = %.3f",
                  s_crino$ratio, s_crino$lo, s_crino$hi,
                  s_crino$p_bf, s_crino$p_perm),
          sprintf("SD ratio %.2f\nBrown-Forsythe P = %.2f",
                  s_macro$ratio, s_macro$p_bf)))

COLG <- c(`Non-diabetic` = "#90a4ae", `Pre-T2D` = "#66bb6a",
          `T2D` = "#1b5e20")

p <- ggplot(long, aes(group, z, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .55, width = .6,
               linewidth = .35) +
  geom_jitter(width = .15, size = .95, alpha = .55, colour = "grey25") +
  geom_text(data = ann, aes(x = 0.55, y = Inf, label = txt),
            inherit.aes = FALSE, hjust = 0, vjust = 1.35,
            size = 2.6, colour = "grey25", lineheight = 1.05) +
  facet_wrap(~ pathway, nrow = 1) +
  scale_fill_manual(values = COLG, guide = "none") +
  coord_cartesian(ylim = c(min(long$z) - .2, max(long$z) + 1.6)) +
  labs(x = NULL,
       y = "Adjusted pathway score (z, within pathway)",
       title = "Dispersion, not location, differs by glycemic status",
       subtitle = sprintf("Discovery cohort, %d donors with complete covariate data",
                          nrow(d2)),
       caption = paste(
         "Scores are residuals after adjustment for donor age, sex, BMI, culture time and islet recovery,",
         "\nstandardised within pathway so the two sets can share an axis; raw score magnitudes differ with",
         "\nset size and are not directly comparable. SD ratios are T2D versus non-diabetic. The crinophagy",
         "\nconfidence interval includes 1.0; the effect is imprecise and sensitive to individual donors",
         "\n(Figure S8C).")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        plot.caption = element_text(size = 6.6, colour = "grey40",
                                    hjust = 0))

ggsave(file.path(FIGS, "Fig1D_dispersion.pdf"), p,
       width = 6.5, height = 4.2, useDingbats = FALSE)

out <- data.frame(
  pathway = c("crinophagy", "comparator"),
  sd_nd = c(s_crino$sd_nd, s_macro$sd_nd),
  sd_t2d = c(s_crino$sd_t2d, s_macro$sd_t2d),
  sd_ratio = c(s_crino$ratio, s_macro$ratio),
  ci_lo = c(s_crino$lo, NA), ci_hi = c(s_crino$hi, NA),
  p_brown_forsythe = c(s_crino$p_bf, s_macro$p_bf),
  p_permutation = c(s_crino$p_perm, NA))
write.csv(out, file.path(FIGS, "Fig1D_values.csv"), row.names = FALSE)

cat("\nWrote figures/Fig1D_dispersion.pdf and Fig1D_values.csv\n")


# =====================================================================
# IF modelling_dataset.csv IS MISSING
#
# It was written by script 02 (02_model_crinophagy.R). To recreate it:
#
#   setwd("/Users/mmani3/Downloads/Human Islets")
#   source("02_model_crinophagy.R")
#
# That script needs only files on disk:
#   proteomics/proc_prot.rds
#   donor.csv
#   frozen_gene_sets_v2_clean.txt
# and the packages GSVA, dplyr, ggplot2, car, broom.
#
# It takes a few minutes (ssGSEA over 6,138 proteins) and writes
# modelling_dataset.csv into the working folder. Then rerun this script.
# =====================================================================
