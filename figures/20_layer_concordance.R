# =====================================================================
# Transcript-protein concordance, discovery cohort
#
# TWO PANELS, TWO QUESTIONS.
#
# Panel A (adapts Wigger et al. Extended Data Fig 4d): do transcript and
#   protein abundance agree across genes? They reported R = 0.3 from
#   5 ND + 5 T2D donors. We have 87 donors with both layers. This is a
#   descriptive context panel: transcript-protein discordance is a known
#   general property, and stating it up front pre-empts the objection
#   that our cross-layer replication should have been tighter.
#
# Panel B (not in their paper, and the one that matters): does the
#   VARIANCE RATIO agree between layers? Our replication rests on
#   UNC13D showing elevated inter-donor variance in proteomics (2.065)
#   and again in transcriptomics (2.001). If variance ratios agree
#   across layers no better than chance genome-wide, that agreement is
#   less impressive than it looks. If they agree modestly and UNC13D
#   sits far out on the diagonal, the cross-platform result is
#   strengthened.
#
#   Either outcome is worth knowing before a reviewer asks.
#
# NOTE: panel A uses the discovery cohort's own RNA-seq (same donors,
# same platform pairing). Panel B compares discovery proteomics against
# discovery transcriptomics - both from HumanIslets, 87 shared donors -
# NOT against GSE164416. That keeps the comparison within one cohort.
# =====================================================================

library(ggplot2)
library(dplyr)
library(patchwork)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

set.seed(1)


# ---- 1. Load both layers --------------------------------------------

prot <- readRDS("./proteomics/proc_prot.rds")
rna  <- readRDS("./Bulk RNA/proc_rnaseq.rds")
donor <- read.csv("donor.csv", stringsAsFactors = FALSE)

cat("Protein matrix:", dim(prot), "\n")
cat("RNA matrix:", dim(rna), "\n")
cat("RNA metadata columns:", paste(colnames(rna)[1:5], collapse = ", "), "\n")

# Proteomics: Entrez-keyed, col 1 is gene_id
p_mat <- as.matrix(prot[, -1, drop = FALSE])
storage.mode(p_mat) <- "numeric"
rownames(p_mat) <- as.character(prot[[1]])

# RNA: first 5 columns are metadata (accession, gene_id, symbol, name,
# biotype). Confirm before slicing.
meta_names <- c("accession","gene_id","symbol","name","biotype")
n_meta <- sum(colnames(rna)[1:6] %in% meta_names)
cat("Metadata columns detected:", n_meta, "\n")

r_meta <- rna[, 1:n_meta]
r_mat  <- as.matrix(rna[, -(1:n_meta), drop = FALSE])
storage.mode(r_mat) <- "numeric"

# Key RNA rows on the same identifier as the proteomics.
# gene_id in proc_rnaseq should be Entrez; verify against overlap.
id_col <- if ("gene_id" %in% colnames(r_meta)) "gene_id" else colnames(r_meta)[2]
rownames(r_mat) <- as.character(r_meta[[id_col]])

shared_genes <- intersect(rownames(p_mat), rownames(r_mat))
cat("Genes in both layers:", length(shared_genes), "\n")

if (length(shared_genes) < 500) {
  cat("\nLow overlap - the two matrices may use different ID systems.\n")
  cat("Protein IDs:", head(rownames(p_mat), 3), "\n")
  cat("RNA IDs    :", head(rownames(r_mat), 3), "\n")
  stop("Resolve identifier mismatch before proceeding.")
}

shared_donors <- intersect(colnames(p_mat), colnames(r_mat))
cat("Donors in both layers:", length(shared_donors), "\n")


# ---- 2. Panel A: abundance concordance ------------------------------

P <- p_mat[shared_genes, shared_donors]
R <- r_mat[shared_genes, shared_donors]

ab <- data.frame(
  gene = shared_genes,
  prot = rowMeans(P, na.rm = TRUE),
  rna  = rowMeans(R, na.rm = TRUE)
) %>% filter(is.finite(prot), is.finite(rna))

r_pearson  <- cor(ab$prot, ab$rna, use = "complete.obs")
r_spearman <- cor(ab$prot, ab$rna, method = "spearman",
                  use = "complete.obs")

cat("\n===== ABUNDANCE CONCORDANCE =====\n")
cat(sprintf("Genes: %d   Donors: %d\n", nrow(ab), length(shared_donors)))
cat(sprintf("Pearson r  = %.3f\n", r_pearson))
cat(sprintf("Spearman r = %.3f\n", r_spearman))
cat("Wigger et al. reported R = 0.3 from 10 samples.\n")

mark_ids <- c(UNC13D = "201294", STX7 = "8417", CTSD = "1509",
              RAB3A = "5864", INS = "3630")
mk_ab <- ab %>% filter(gene %in% mark_ids) %>%
  mutate(label = names(mark_ids)[match(gene, mark_ids)])

pa <- ggplot(ab, aes(prot, rna)) +
  geom_point(alpha = .18, size = .7, colour = "#4a90d9") +
  geom_smooth(method = "lm", colour = "grey25", linewidth = .5,
              se = FALSE) +
  geom_point(data = mk_ab, colour = "#c62828", size = 2) +
  ggrepel::geom_text_repel(data = mk_ab, aes(label = label),
                           size = 2.7, seed = 1, box.padding = .4) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.6,
           label = sprintf("r = %.2f", r_pearson), size = 3.2) +
  labs(x = "Protein abundance (mean, processed)",
       y = "Transcript abundance (mean, processed)",
       title = "A   Abundance concordance across layers",
       subtitle = sprintf("%s genes, %d donors with both layers",
                          format(nrow(ab), big.mark = ","),
                          length(shared_donors))) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"))


# ---- 3. Panel B: variance-ratio concordance -------------------------
# The question our replication actually depends on.

md <- donor[match(shared_donors, donor$record_id), ]
rownames(md) <- md$record_id
md$group <- factor(md$diagnosis_computed,
                   levels = c("None", "Pre.T2D", "Type2"))

COV <- c("donorage","donorsex","bodymassindex",
         "predistributionculturetime","percentieqrecoverypostculture")
ok <- complete.cases(md[, COV]) & md$group %in% c("None","Type2")
md2 <- droplevels(md[ok, ])
md2$donorsex <- factor(md2$donorsex)

cat("\n===== VARIANCE-RATIO CONCORDANCE =====\n")
cat("Donors used (ND vs T2D, complete covariates):", nrow(md2), "\n")
print(table(md2$group))

if (nrow(md2) < 30) stop("Too few donors for a stable comparison.")

form <- as.formula(paste("y ~", paste(COV, collapse = " + ")))

vratio <- function(mat) {
  ids <- rownames(md2)
  apply(mat[, ids, drop = FALSE], 1, function(y) {
    if (!is.finite(sd(y)) || sd(y) == 0) return(NA_real_)
    r <- residuals(lm(form, data = transform(md2, y = y)))
    s1 <- sd(r[md2$group == "Type2"]); s0 <- sd(r[md2$group == "None"])
    if (!is.finite(s0) || s0 == 0) return(NA_real_)
    s1 / s0
  })
}

cat("Computing variance ratios in both layers...\n")
vr_p <- vratio(P)
vr_r <- vratio(R)

vc <- data.frame(gene = shared_genes, prot = vr_p, rna = vr_r) %>%
  filter(is.finite(prot), is.finite(rna),
         prot > 0.1, prot < 10, rna > 0.1, rna < 10)

cv_p <- cor(log2(vc$prot), log2(vc$rna))
cv_s <- cor(vc$prot, vc$rna, method = "spearman")

cat(sprintf("Genes compared: %d\n", nrow(vc)))
cat(sprintf("Pearson r (log2 ratios) = %.3f\n", cv_p))
cat(sprintf("Spearman r              = %.3f\n", cv_s))

mk_vc <- vc %>% filter(gene %in% mark_ids) %>%
  mutate(label = names(mark_ids)[match(gene, mark_ids)])
cat("\nMarked genes:\n")
print(mk_vc %>% mutate(across(c(prot, rna), ~round(.x, 3))),
      row.names = FALSE)

cat("\nHOW TO READ THIS:\n")
if (cv_s > 0.15) {
  cat("  Variance ratios agree modestly across layers. UNC13D appearing\n")
  cat("  high in both is therefore not merely expected; check where it\n")
  cat("  sits relative to the bulk of genes.\n")
} else {
  cat("  Variance ratios agree poorly across layers genome-wide. Our\n")
  cat("  cross-platform UNC13D result is therefore NOT explained by a\n")
  cat("  general tendency for the two layers to agree, which makes it\n")
  cat("  more striking - but it also means cross-layer agreement is a\n")
  cat("  weak criterion in general. State this plainly.\n")
}

pb <- ggplot(vc, aes(log2(prot), log2(rna))) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(alpha = .15, size = .7, colour = "#4a90d9") +
  geom_point(data = mk_vc, colour = "#c62828", size = 2) +
  ggrepel::geom_text_repel(data = mk_vc, aes(label = label),
                           size = 2.7, seed = 1, box.padding = .4) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.2, vjust = 1.6,
           label = sprintf("r = %.2f", cv_p), size = 3.2) +
  labs(x = expression(paste(log[2], " variance ratio, protein")),
       y = expression(paste(log[2], " variance ratio, transcript")),
       title = "B   Variance-ratio concordance across layers",
       subtitle = sprintf("Same %d donors, T2D vs non-diabetic, covariate adjusted",
                          nrow(md2)),
       caption = paste(
         "Panel A adapts the transcript-protein comparison of Wigger et al. (Nat Metab 2021),",
         "\nwho reported R = 0.3 from ten samples. Panel B asks the same question of variance",
         "\nratios rather than abundances; to our knowledge this has not been reported.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        plot.caption = element_text(size = 6.8, colour = "grey40", hjust = 0))

fig <- pa | pb
ggsave(file.path(FIGS, "FigS_layer_concordance.pdf"), fig,
       width = 9, height = 4.5)

write.csv(vc, file.path(FIGS, "layer_concordance.csv"), row.names = FALSE)

cat("\nWrote FigS_layer_concordance.pdf and layer_concordance.csv\n")
