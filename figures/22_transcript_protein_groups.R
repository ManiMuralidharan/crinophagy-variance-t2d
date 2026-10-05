# =====================================================================
# Transcript vs protein, panel genes, all glycemic groups
#
# WHAT THIS SHOWS. For every gene in the frozen sets, mean transcript
# abundance against mean protein abundance, computed separately in
# non-diabetic, pre-T2D and T2D donors. Adapts the layer-comparison
# format of Wigger et al. (Nat Metab 2021, Extended Data Fig 4d), who
# plotted all genes from 10 samples; here restricted to the pre-registered
# panel and split by glycemic status, in the 87 discovery donors with
# both layers.
#
# WHAT IT IS FOR. Two things:
#   1. Context. A gene sitting far off the diagonal is measured very
#      differently by the two platforms, which matters when a result
#      replicates across layers (UNC13D) or fails to (RAB3A).
#   2. Stability. If a gene's position shifts between groups, its
#      transcript-protein relationship is itself disease-dependent,
#      which is worth knowing before interpreting either layer alone.
#
# WHAT IT IS NOT. This plots abundance, not variance. It does not test
# anything and no statistic from it enters the manuscript claims. It is
# a descriptive supplementary panel.
# =====================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)


# ---- 1. Load ---------------------------------------------------------

prot  <- readRDS("./proteomics/proc_prot.rds")
rna   <- readRDS("./Bulk RNA/proc_rnaseq.rds")
donor <- read.csv("donor.csv", stringsAsFactors = FALSE)
sets  <- read.delim("frozen_gene_sets_v2_clean.txt",
                    stringsAsFactors = FALSE)

p_mat <- as.matrix(prot[, -1, drop = FALSE])
storage.mode(p_mat) <- "numeric"
rownames(p_mat) <- as.character(prot[[1]])

meta_names <- c("accession","gene_id","symbol","name","biotype")
n_meta <- sum(colnames(rna)[1:6] %in% meta_names)
r_meta <- rna[, 1:n_meta]
r_mat  <- as.matrix(rna[, -(1:n_meta), drop = FALSE])
storage.mode(r_mat) <- "numeric"
id_col <- if ("gene_id" %in% colnames(r_meta)) "gene_id" else colnames(r_meta)[2]
rownames(r_mat) <- as.character(r_meta[[id_col]])

shared <- intersect(colnames(p_mat), colnames(r_mat))
cat("Donors with both layers:", length(shared), "\n")

md <- donor[match(shared, donor$record_id), ]
rownames(md) <- md$record_id
md$group <- factor(md$diagnosis_computed,
                   levels = c("None", "Pre.T2D", "Type2"),
                   labels = c("Non-diabetic", "Pre-T2D", "T2D"))
print(table(md$group))


# ---- 2. Panel genes, plus post-hoc comparators -----------------------

key <- sets %>%
  dplyr::select(SYMBOL, ENTREZID, set) %>%
  mutate(ENTREZID = as.character(ENTREZID))

extra <- data.frame(SYMBOL = c("STX1A","UNC13A","UNC13B"),
                    ENTREZID = c("6804","23025","10497"),
                    set = "post-hoc", stringsAsFactors = FALSE)
key <- bind_rows(key, extra) %>% distinct(SYMBOL, .keep_all = TRUE)

key <- key %>%
  filter(ENTREZID %in% rownames(p_mat), ENTREZID %in% rownames(r_mat))

cat("Panel genes measured in both layers:", nrow(key), "of", nrow(sets) + 3, "\n")
missing <- setdiff(c(sets$SYMBOL, extra$SYMBOL), key$SYMBOL)
if (length(missing))
  cat("Not in both layers:", paste(missing, collapse = ", "), "\n")


# ---- 3. Per-group means ----------------------------------------------

grp_mean <- function(mat, ids, donors) {
  sapply(levels(md$group), function(g) {
    d <- donors[md$group[match(donors, rownames(md))] == g]
    rowMeans(mat[ids, d, drop = FALSE], na.rm = TRUE)
  })
}

pm <- grp_mean(p_mat, key$ENTREZID, shared)
rm_ <- grp_mean(r_mat, key$ENTREZID, shared)

df <- bind_rows(lapply(colnames(pm), function(g) {
  data.frame(SYMBOL = key$SYMBOL, set = key$set, group = g,
             prot = pm[, g], rna = rm_[, g],
             stringsAsFactors = FALSE)
})) %>%
  filter(is.finite(prot), is.finite(rna), prot > 0, rna > 0) %>%
  mutate(group = factor(group, levels = levels(md$group)),
         lp = log10(prot), lr = log10(rna))

cat("\nPoints plotted:", nrow(df), "\n")


# ---- 4. Correlation within each group --------------------------------

cors <- df %>%
  group_by(group) %>%
  summarise(n = n(),
            pearson = cor(lp, lr),
            spearman = cor(prot, rna, method = "spearman"),
            .groups = "drop")

cat("\n===== TRANSCRIPT-PROTEIN CORRELATION, PANEL GENES =====\n")
print(as.data.frame(cors), digits = 3, row.names = FALSE)
cat("\nWigger et al. reported R = 0.3 across all genes in 10 samples.\n")

# Which genes sit furthest off the common trend?
fit <- lm(lr ~ lp, data = df)
df$resid <- residuals(fit)
cat("\nFurthest from the transcript-protein trend (T2D group):\n")
print(df %>% filter(group == "T2D") %>%
        arrange(desc(abs(resid))) %>%
        dplyr::select(SYMBOL, set, lp, lr, resid) %>%
        head(8), digits = 3, row.names = FALSE)


# ---- 5. Does a gene's position shift between groups? -----------------

shift <- df %>%
  dplyr::select(SYMBOL, group, lp, lr) %>%
  pivot_wider(names_from = group, values_from = c(lp, lr)) %>%
  mutate(d_prot = `lp_T2D` - `lp_Non-diabetic`,
         d_rna  = `lr_T2D` - `lr_Non-diabetic`) %>%
  arrange(desc(abs(d_prot) + abs(d_rna)))

cat("\nLargest T2D vs non-diabetic shifts in either layer:\n")
print(shift %>% dplyr::select(SYMBOL, d_prot, d_rna) %>% head(8),
      digits = 3, row.names = FALSE)


# ---- 6. Figure -------------------------------------------------------

lab_genes <- c("UNC13D","STX7","CTSD","RAB3A","RAB27A","INS",
               "CHGA","LAMP1","STX1A","UNC13B")

fig <- ggplot(df, aes(lp, lr)) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey70",
              linewidth = .4) +
  geom_point(aes(colour = set), size = 1.9, alpha = .85) +
  geom_text_repel(data = df %>% filter(SYMBOL %in% lab_genes),
                  aes(label = SYMBOL), size = 2.4, seed = 1,
                  box.padding = .35, max.overlaps = 25,
                  segment.size = .2, colour = "grey20") +
  facet_wrap(~ group, nrow = 1) +
  scale_colour_manual(values = c(crinophagy = "#1b5e20",
                                 macroautophagy = "#90a4ae",
                                 `post-hoc` = "#ef6c00"),
                      name = NULL) +
  labs(x = expression(paste(log[10], " protein abundance")),
       y = expression(paste(log[10], " transcript abundance")),
       title = "Transcript and protein abundance for the pre-registered panel",
       subtitle = sprintf("%d donors with both layers (%s)",
                          length(shared),
                          paste(sprintf("%s n=%d", levels(md$group),
                                        as.integer(table(md$group))),
                                collapse = ", ")),
       caption = paste(
         "Group means per gene. Grey line: linear fit across all points. Descriptive only;",
         "\nno statistic from this panel enters the manuscript claims. Format adapted from",
         "\nWigger et al. (Nat Metab 2021); data and analysis are ours.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(size = 9, face = "bold"),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        plot.caption = element_text(size = 6.8, colour = "grey40", hjust = 0),
        legend.position = "top")

ggsave(file.path(FIGS, "FigS_transcript_protein_bygroup.pdf"), fig,
       width = 10, height = 4.5)

write.csv(df, file.path(FIGS, "transcript_protein_bygroup.csv"),
          row.names = FALSE)

cat("\nWrote FigS_transcript_protein_bygroup.pdf and the underlying table\n")
