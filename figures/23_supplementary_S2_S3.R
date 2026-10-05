# =====================================================================
# Supplementary S2 and S3
#
#   S2  Permutation null distributions behind the set-level tests.
#       The p-values are already reported; this shows the nulls they
#       came from, so a reader can see the observed statistic against
#       its reference rather than taking the number on trust.
#
#   S3  Donor-level view of the discovery cohort, including the
#       leave-one-out sensitivity. This is the figure that shows the
#       pathway-level result is fragile. We include it deliberately:
#       the fragility is stated in the Results and the Discussion, and
#       a reviewer who reads that should be able to see it.
#
# Requires for S2: mat, pd2, bg1, bg2 (scripts 09 + 13)
# Requires for S3: discovery proteomics and donor.csv
# =====================================================================

library(ggplot2)
library(dplyr)
library(patchwork)
library(car)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

set.seed(1)
NPERM <- 10000


# =====================================================================
# S2  Permutation distributions
# =====================================================================

stopifnot(exists("bg1"), exists("bg2"))

sets <- read.delim("frozen_gene_sets_v2_clean.txt",
                   stringsAsFactors = FALSE)

# variance ratios for set members, per model, from the frozen table
full <- read.csv(file.path("replication_GSE164416", "canonical_full.csv"),
                 stringsAsFactors = FALSE)

perm_panel <- function(model_tag, bg, title_prefix) {
  res <- full %>% dplyr::filter(model == model_tag)

  out <- lapply(c("crinophagy", "comparator"), function(s) {
    obs_g <- res$sd_ratio[res$set == s]
    obs_g <- obs_g[is.finite(obs_g)]
    n <- length(obs_g)
    obs <- median(obs_g)
    null <- replicate(NPERM, median(sample(bg, n)))
    p <- (sum(null >= obs) + 1) / (NPERM + 1)
    list(df = data.frame(set = s, null = null),
         obs = data.frame(set = s, obs = obs, n = n, p = p))
  })

  nulls <- bind_rows(lapply(out, `[[`, "df"))
  obss  <- bind_rows(lapply(out, `[[`, "obs"))

  nulls$set <- factor(nulls$set, levels = c("crinophagy", "comparator"),
                      labels = c("Crinophagy (pre-registered)",
                                 "Macroautophagy comparator"))
  obss$set <- factor(obss$set, levels = c("crinophagy", "comparator"),
                     labels = levels(nulls$set))
  obss$lab <- sprintf("observed %.3f\nP = %.3f", obss$obs, obss$p)

  print(obss)

  ggplot(nulls, aes(null)) +
    geom_histogram(bins = 60, fill = "grey85", colour = NA) +
    geom_vline(data = obss, aes(xintercept = obs), colour = "#1b5e20",
               linewidth = .7) +
    geom_text(data = obss, aes(x = obs, y = Inf, label = lab),
              hjust = -0.08, vjust = 1.4, size = 2.6, colour = "#1b5e20") +
    facet_wrap(~ set, nrow = 1, scales = "free_y") +
    labs(x = "Median variance ratio of a random gene set of equal size",
         y = "Permutations",
         title = paste0(title_prefix,
                        "   Set-level permutation null"),
         subtitle = sprintf("%s random sets drawn from the genome-wide background",
                            format(NPERM, big.mark = ","))) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 8.5, face = "bold"),
          plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8, colour = "grey35"))
}

cat("===== S2: MODEL 1 (batch only) =====\n")
s2a <- perm_panel("M1_batch", bg1, "A")

cat("\n===== S2: MODEL 2 (fully adjusted) =====\n")
s2b <- perm_panel("M2_full", bg2, "B")

s2 <- s2a / s2b +
  plot_annotation(
    caption = paste(
      "Both gene sets were frozen with a timestamped DOI before analysis, so these tests are",
      "\npre-specified. The set-level evidence is stronger under batch-only adjustment; the fully",
      "\nadjusted model (B) is primary because the confounds it addresses are present in these data.")
  ) & theme(plot.caption = element_text(size = 6.8, colour = "grey40",
                                        hjust = 0))

ggsave(file.path(FIGS, "FigS2_permutation_nulls.pdf"), s2,
       width = 8, height = 6.5)
cat("\nWrote FigS2_permutation_nulls.pdf\n")


# =====================================================================
# S3  Discovery cohort, donor level, and leave-one-out
# =====================================================================

prot  <- readRDS("./proteomics/proc_prot.rds")
donor <- read.csv("donor.csv", stringsAsFactors = FALSE)

# reuse the frozen modelling dataset if present, else stop with a note
if (file.exists("modelling_dataset.csv")) {
  d <- read.csv("modelling_dataset.csv", na.strings = c("NA",""),
                stringsAsFactors = FALSE)
  # "None" must not be read as missing
  d$diagnosis_computed[is.na(d$diagnosis_computed)] <- "None"
} else {
  stop("modelling_dataset.csv not found - rerun script 02 first.")
}

d$group <- factor(d$diagnosis_computed,
                  levels = c("None","Pre.T2D","Type2"),
                  labels = c("Non-diabetic","Pre-T2D","T2D"))

COV <- c("donorage","donorsex","bodymassindex",
         "predistributionculturetime","percentieqrecoverypostculture")
ok <- complete.cases(d[, c("crino_raw","macro_raw", COV)])
d2 <- droplevels(d[ok, ])
d2$donorsex <- factor(d2$donorsex)

form <- as.formula(paste("crino_raw ~", paste(COV, collapse = " + ")))
d2$resid <- residuals(lm(form, data = d2))
d2$resid_mac <- residuals(lm(as.formula(paste("macro_raw ~",
                        paste(COV, collapse = " + "))), data = d2))

cat("\n===== S3: DISCOVERY COHORT =====\n")
print(d2 %>% group_by(group) %>%
        summarise(n = n(), sd = sd(resid), .groups = "drop") %>%
        as.data.frame(), digits = 4, row.names = FALSE)

lev <- leveneTest(resid ~ group, data = d2)
cat(sprintf("Brown-Forsythe P = %.4f\n", lev$`Pr(>F)`[1]))

COLD <- c(`Non-diabetic` = "#90a4ae", `Pre-T2D` = "#66bb6a",
          `T2D` = "#1b5e20")

# --- A: donor-level adjusted scores, both pathways ---
long <- bind_rows(
  d2 %>% transmute(group, score = resid, pathway = "Crinophagy"),
  d2 %>% transmute(group, score = resid_mac, pathway = "Macroautophagy")
) %>%
  group_by(pathway) %>%
  mutate(z = (score - mean(score)) / sd(score)) %>%
  ungroup()

pa <- ggplot(long, aes(group, z, fill = group)) +
  geom_boxplot(outlier.shape = NA, alpha = .55, width = .6,
               linewidth = .3) +
  geom_jitter(width = .16, size = .9, alpha = .55, colour = "grey25") +
  facet_wrap(~ pathway, nrow = 1) +
  scale_fill_manual(values = COLD, guide = "none") +
  labs(x = NULL, y = "Adjusted score (z, within pathway)",
       title = "A   Donor-level pathway scores, discovery cohort") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 8.5, face = "bold"),
        plot.title = element_text(size = 10, face = "bold"))

# --- B: score against HbA1c ---
pb <- ggplot(d2 %>% filter(!is.na(hba1c)), aes(hba1c, resid)) +
  geom_point(aes(colour = group), size = 1.6, alpha = .85) +
  geom_smooth(method = "lm", colour = "grey30", linewidth = .5,
              se = TRUE, fill = "grey88") +
  scale_colour_manual(values = COLD, name = NULL) +
  labs(x = "HbA1c (%)", y = "Adjusted crinophagy score",
       title = "B   No mean association with glycemia",
       subtitle = "beta = 11.8, 95% CI -1.7 to 25.2, P = 0.086") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        legend.position = "top")

# --- C: leave-one-out over T2D donors ---
t2 <- which(d2$group == "T2D")
loo <- do.call(rbind, lapply(t2, function(i) {
  s <- d2[-i, ]
  data.frame(dropped = d2$record_id[i],
             p = leveneTest(resid ~ group, data = droplevels(s))$`Pr(>F)`[1])
})) %>% arrange(desc(p))
loo$dropped <- factor(loo$dropped, levels = rev(loo$dropped))

cat("\nLeave-one-out worst case: P =", round(max(loo$p), 4), "\n")

pc <- ggplot(loo, aes(p, dropped)) +
  geom_vline(xintercept = 0.05, linetype = "dashed", colour = "#c62828") +
  geom_point(size = 1.9, colour = "#1b5e20") +
  labs(x = "Brown-Forsythe P with that donor removed", y = NULL,
       title = "C   Leave-one-out over T2D donors",
       subtitle = sprintf("Worst case P = %.3f; %d of %d exclusions cross 0.05",
                          max(loo$p), sum(loo$p > 0.05), nrow(loo)),
       caption = paste(
         "The discovery pathway-level result is sensitive to individual donors. This is why the",
         "\nreplication cohort and the genome-wide background, rather than this test alone, carry",
         "\nthe claim. Shown here so the fragility stated in the text can be seen directly.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 7),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        plot.caption = element_text(size = 6.8, colour = "grey40", hjust = 0))

s3 <- (pa | pb) / pc + plot_layout(heights = c(1, 1.1))

ggsave(file.path(FIGS, "FigS3_discovery_donor_level.pdf"), s3,
       width = 9, height = 8)

write.csv(loo, file.path(FIGS, "leave_one_out.csv"), row.names = FALSE)

cat("\nWrote FigS3_discovery_donor_level.pdf and leave_one_out.csv\n")
