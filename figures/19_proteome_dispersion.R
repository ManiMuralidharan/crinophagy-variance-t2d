# =====================================================================
# Whole-proteome dispersion, discovery cohort
#
# WHAT THIS ADAPTS. Wigger et al. (Nat Metab 2021) observed that islet
# proteomes from ND donors clustered tightly in PCA while T2D proteomes
# scattered. That observation rests on 5 ND and 5 T2D samples and was
# reported qualitatively, with no test.
#
# WHAT THIS DOES DIFFERENTLY. Same question, 129 donors, covariate
# adjusted, and tested:
#   - PCA on the adjusted proteome
#   - mean Euclidean distance from each group's centroid in PC space
#   - permutation test on the ND vs T2D difference in that distance
#   - sensitivity to the number of PCs retained
#
# This converts a qualitative observation into a quantitative result in
# a cohort 13x larger, and it gives the paper a whole-proteome panel to
# sit before the gene-level results.
#
# NOTE ON FRAMING. This is an adaptation of their analytical approach,
# not a reuse of their figure. Cite them for the observation; the data
# and the test are ours.
#
# Output: two panels intended for Figure 1.
# =====================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

set.seed(1)
NPERM <- 10000

GROUPS <- c("NoDiabetes", "Pre.T2D", "Type2")
LABS   <- c(NoDiabetes = "Non-diabetic", Pre.T2D = "Pre-T2D",
            Type2 = "T2D")
COL    <- c(`Non-diabetic` = "#90a4ae", `Pre-T2D` = "#66bb6a",
            `T2D` = "#1b5e20")


# ---- 1. Load discovery data -----------------------------------------

prot <- readRDS("./proteomics/proc_prot.rds")
donor <- read.csv("donor.csv", stringsAsFactors = FALSE)

expr <- as.matrix(prot[, -1, drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- as.character(prot[[1]])

md <- donor[match(colnames(expr), donor$record_id), ]
rownames(md) <- md$record_id
md$group <- factor(md$diagnosis_computed,
                   levels = c("None", "Pre.T2D", "Type2"),
                   labels = GROUPS)

cat("Donors:", ncol(expr), " Proteins:", nrow(expr), "\n")
print(table(md$group))


# ---- 2. Covariate adjustment ----------------------------------------
# Their analysis did not adjust. Ours must: if T2D donors differ in age,
# BMI or culture time, unadjusted scatter could reflect those rather
# than disease. Adjusting first makes the dispersion claim meaningful.

COV <- c("donorage", "donorsex", "bodymassindex",
         "predistributionculturetime", "percentieqrecoverypostculture")

keep <- complete.cases(md[, COV]) & !is.na(md$group)
md2  <- droplevels(md[keep, ])
e2   <- expr[, rownames(md2)]

cat("\nAfter dropping donors with missing covariates:", ncol(e2), "\n")
print(table(md2$group))

md2$donorsex <- factor(md2$donorsex)

# Residualise every protein on the covariates (not on group)
form <- as.formula(paste("y ~", paste(COV, collapse = " + ")))
resid_mat <- t(apply(e2, 1, function(y) {
  residuals(lm(form, data = transform(md2, y = y)))
}))
colnames(resid_mat) <- colnames(e2)

cat("Adjusted matrix:", dim(resid_mat), "\n")


# ---- 3. PCA ---------------------------------------------------------

pca <- prcomp(t(resid_mat), center = TRUE, scale. = FALSE)
ve <- 100 * pca$sdev^2 / sum(pca$sdev^2)

cat("\nVariance explained by first 5 PCs:\n")
print(round(ve[1:5], 1))

sc <- as.data.frame(pca$x[, 1:10])
sc$group <- md2$group
sc$label <- LABS[as.character(sc$group)]
sc$label <- factor(sc$label, levels = LABS)


# ---- 4. Centroid distance, with a test ------------------------------
# Their claim was visual. This is the quantity behind it.

centroid_dist <- function(scores, grp, npc) {
  X <- as.matrix(scores[, paste0("PC", 1:npc), drop = FALSE])
  out <- numeric(nrow(X))
  for (g in levels(grp)) {
    i <- which(grp == g)
    if (length(i) < 2) next
    cen <- colMeans(X[i, , drop = FALSE])
    out[i] <- sqrt(rowSums(sweep(X[i, , drop = FALSE], 2, cen)^2))
  }
  out
}

NPC <- 5
sc$dist <- centroid_dist(sc, sc$group, NPC)

cat("\n===== MEAN DISTANCE FROM GROUP CENTROID (first", NPC, "PCs) =====\n")
summ <- sc %>% group_by(label) %>%
  summarise(n = n(), mean_dist = mean(dist), median_dist = median(dist),
            .groups = "drop")
print(as.data.frame(summ), digits = 4, row.names = FALSE)

# Permutation test: ND vs T2D, labels shuffled between those two groups
sub <- sc %>% filter(group %in% c("NoDiabetes", "Type2")) %>%
  droplevels()

obs <- mean(sub$dist[sub$group == "Type2"]) -
       mean(sub$dist[sub$group == "NoDiabetes"])

null <- numeric(NPERM)
g0 <- sub$group
for (i in seq_len(NPERM)) {
  gp <- sample(g0)
  # centroids must be recomputed under the permuted labels
  d <- centroid_dist(sub, gp, NPC)
  null[i] <- mean(d[gp == "Type2"]) - mean(d[gp == "NoDiabetes"])
}

p_perm <- (sum(abs(null) >= abs(obs)) + 1) / (NPERM + 1)

cat(sprintf("\nObserved difference (T2D - ND): %.2f\n", obs))
cat(sprintf("Permutation p (%s shuffles): %.4f\n",
            format(NPERM, big.mark = ","), p_perm))
cat(sprintf("Ratio of mean distances: %.3f\n",
            mean(sub$dist[sub$group == "Type2"]) /
            mean(sub$dist[sub$group == "NoDiabetes"])))

# Also the standard multivariate dispersion test, if available
if (requireNamespace("vegan", quietly = TRUE)) {
  library(vegan)
  D <- dist(as.matrix(sc[, paste0("PC", 1:NPC)]))
  bd <- betadisper(D, sc$group)
  cat("\nPERMANOVA dispersion test (vegan::betadisper):\n")
  print(permutest(bd, permutations = 999))
} else {
  cat("\nInstall vegan for the betadisper cross-check.\n")
}


# ---- 5. Sensitivity to the number of PCs ----------------------------
# A result that only holds at one arbitrary choice of NPC is not a
# result. This reports the ratio across a range.

cat("\n===== SENSITIVITY TO PC COUNT =====\n")
for (k in c(2, 3, 5, 10)) {
  d <- centroid_dist(sc, sc$group, k)
  r <- mean(d[sc$group == "Type2"]) / mean(d[sc$group == "NoDiabetes"])
  cat(sprintf("  %2d PCs: T2D/ND mean distance ratio = %.3f\n", k, r))
}


# ---- 6. Panels ------------------------------------------------------

pa <- ggplot(sc, aes(PC1, PC2, colour = label)) +
  stat_ellipse(aes(fill = label), geom = "polygon", alpha = .12,
               colour = NA, level = .68) +
  geom_point(size = 1.8, alpha = .85) +
  scale_colour_manual(values = COL, name = NULL) +
  scale_fill_manual(values = COL, guide = "none") +
  labs(x = sprintf("PC1 (%.1f%%)", ve[1]),
       y = sprintf("PC2 (%.1f%%)", ve[2]),
       title = "A   Islet proteome, covariate-adjusted",
       subtitle = sprintf("%d donors. Ellipses: 68%% contour.", nrow(sc))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        legend.position = "top")

pb <- ggplot(sc, aes(label, dist, fill = label)) +
  geom_boxplot(outlier.shape = NA, alpha = .55, width = .55,
               linewidth = .3) +
  geom_jitter(width = .16, size = 1, alpha = .6, colour = "grey25") +
  scale_fill_manual(values = COL, guide = "none") +
  labs(x = NULL, y = sprintf("Distance from group centroid (%d PCs)", NPC),
       title = "B   Proteome dispersion by glycemic status",
       subtitle = sprintf("T2D/ND ratio %.2f, permutation P = %.3f",
                          mean(sub$dist[sub$group == "Type2"]) /
                          mean(sub$dist[sub$group == "NoDiabetes"]),
                          p_perm),
       caption = paste(
         "Adjusted for age, sex, BMI, culture time and islet recovery before PCA.",
         "\nApproach adapted from Wigger et al. (Nat Metab 2021), who reported this",
         "\nqualitatively in 5 ND and 5 T2D proteomes; here tested in the full cohort.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        plot.caption = element_text(size = 6.8, colour = "grey40", hjust = 0))

fig <- pa | pb

ggsave(file.path(FIGS, "Fig1_proteome_dispersion.pdf"), fig,
       width = 9, height = 4.5)

write.csv(sc[, c("PC1","PC2","PC3","group","dist")],
          file.path(FIGS, "proteome_dispersion.csv"))

cat("\nWrote Fig1_proteome_dispersion.pdf and proteome_dispersion.csv\n")
