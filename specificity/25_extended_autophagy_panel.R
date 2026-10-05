# =====================================================================
# EXPLORATORY: extended autophagy / lysosomal fusion panel
#
# STATUS, STATED BEFORE ANYTHING ELSE
#
#   This panel is POST-HOC. It was assembled after the primary analysis
#   and after the GABARAP result was known. It is NOT part of the
#   pre-registered design, it does NOT replace the frozen 14-gene
#   comparator, and no set-level p-value from it should be presented
#   alongside the pre-registered tests as though equivalent.
#
#   The frozen sets and their results stand exactly as reported:
#     crinophagy   31 genes, median 1.178, permutation P = 0.036
#     comparator   14 genes, median 1.036, permutation P = 0.72
#
# THE QUESTION THIS ASKS
#
#   Where in the autophagy-lysosome system does elevated inter-donor
#   variance sit? Specifically, is it confined to components acting at
#   granule-lysosome fusion, or does it extend upstream to initiation
#   and conjugation machinery?
#
# THE OVERLAP PROBLEM, HANDLED EXPLICITLY
#
#   Several genes proposed for a "shared fusion machinery" group are
#   already members of the frozen crinophagy set: STX7, VAMP7, VAMP8,
#   VPS11, VPS16, VPS33A, PIKFYVE. STX7 is the second-ranked result in
#   the primary analysis. A comparator sharing members with the test
#   set cannot function as a comparator. Those genes are therefore
#   flagged as "in frozen crinophagy set" and are EXCLUDED from every
#   group-level summary computed here, while still being plotted so the
#   reader can see where they fall.
#
# Requires: mat, pd2, bg2 in session.
# =====================================================================

library(ggplot2)
library(dplyr)
library(car)
library(org.Hs.eg.db)
library(AnnotationDbi)

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

stopifnot(exists("mat"), exists("pd2"), exists("bg2"))
set.seed(1)


# ---- 1. The extended panel, grouped by pathway position -------------

panel <- list(
  `Upstream signalling` = c(
    "MTOR","RPTOR","PRKAA1","PRKAA2","TBK1","STK4","IKBKB","TMEM173",
    "SIRT3","SIRT4","SIRT5","SIRT6","KAT2B","CREBBP","HDAC2"),

  `Initiation (ULK / PI3KC3)` = c(
    "ULK1","ULK2","ATG13","ATG101","RB1CC1",
    "PIK3C3","BECN1","ATG14","UVRAG","NRBF2","RUBCN","RUBCNL",
    "TRIM22"),

  `Early membrane / omegasome` = c(
    "ATG9A","ATG2A","ATG2B","TMEM49","TMEM41B","VPS13A","VPS13C",
    "WIPI2","ZFYVE1"),

  `Conjugation` = c(
    "ATG3","ATG5","ATG7","ATG10","ATG12","ATG16L1","ATG4B"),

  `ATG8 family` = c(
    "MAP1LC3A","MAP1LC3B","GABARAP","GABARAPL1","GABARAPL2"),

  `Cargo receptors` = c(
    "SQSTM1","NBR1","CALCOCO2","BNIP3L","TAX1BP1"),

  `ESCRT / closure` = c(
    "VPS4A","VPS4B","VPS37A","CHMP2A","VPS25","CHMP4B","VPS28",
    "TSG101"),

  `Rab / effectors (fusion)` = c(
    "RAB7A","RAB2A","RAB39A","RAB21","RAB24","RAB33B","RAB32","RAB38",
    "ARL8A","ARL8B","TBC1D25","RILP","FYCO1","OSBPL1A"),

  `Tethering (fusion)` = c(
    "PLEKHM1","VPS11","VPS16","VPS18","VPS33A","VPS39","VPS41",
    "EPG5","TECPR1","BIRC6"),

  `SNAREs (fusion)` = c(
    "VAMP8","VAMP7","STX17","SNAP29","SNAP47","YKT6","STX7","VTI1B",
    "SCFD1"),

  `Lipid enzymes (fusion)` = c(
    "PIKFYVE","INPP5E","PI4K2A","PIP5K1C"),

  `Amphisome / endosome` = c(
    "VPS4A","RAB11A","VAMP3","SNAP23","ZFYVE26","SNAPIN",
    "GORASP2","NPC1","C9orf72","CCZ1","MON1A")
)

pan <- bind_rows(lapply(names(panel), function(g)
  data.frame(SYMBOL = panel[[g]], group = g, stringsAsFactors = FALSE)
)) %>% distinct(SYMBOL, .keep_all = TRUE)

cat("Genes proposed:", nrow(pan), "\n")


# ---- 2. Flag overlap with the frozen sets ---------------------------

frozen <- read.delim("frozen_gene_sets_v2_clean.txt",
                     stringsAsFactors = FALSE)
crino_frozen <- frozen$SYMBOL[frozen$set == "crinophagy"]
comp_frozen  <- frozen$SYMBOL[frozen$set == "comparator"]

pan <- pan %>%
  mutate(status = case_when(
    SYMBOL %in% crino_frozen ~ "In frozen crinophagy set",
    SYMBOL %in% comp_frozen  ~ "In frozen comparator set",
    TRUE                     ~ "New (post-hoc)"))

cat("\n===== OVERLAP WITH FROZEN SETS =====\n")
print(table(pan$status))
cat("\nOverlapping with crinophagy set:\n  ",
    paste(pan$SYMBOL[pan$status == "In frozen crinophagy set"],
          collapse = ", "), "\n")
cat("\nThese are plotted but EXCLUDED from all group summaries below.\n")


# ---- 3. Map and test ------------------------------------------------

mp <- AnnotationDbi::select(org.Hs.eg.db, keys = pan$SYMBOL,
                            keytype = "SYMBOL", columns = "ENSEMBL")
mp <- mp[!is.na(mp$ENSEMBL) & mp$ENSEMBL %in% rownames(mat), ]
mp <- mp[!duplicated(mp$SYMBOL), ]

pan <- pan %>% inner_join(mp, by = "SYMBOL") %>% rename(ensembl = ENSEMBL)

cat("\nMeasured in the replication matrix:", nrow(pan), "\n")
absent <- setdiff(unlist(panel), pan$SYMBOL)
cat("Not measured:", length(absent), "\n  ",
    paste(head(absent, 25), collapse = ", "),
    if (length(absent) > 25) "..." else "", "\n")

test_one <- function(ens) {
  d <- pd2
  d$y <- mat[ens, rownames(d)]
  d$r <- residuals(lm(y ~ batch + imm + beta, data = d))
  a <- d$r[d$group == "T2D"]; b <- d$r[d$group == "ND"]
  boot <- replicate(2000,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))
  c(ratio = sd(a)/sd(b),
    lo = unname(quantile(boot, .025)), hi = unname(quantile(boot, .975)),
    p = leveneTest(r ~ group, data = d)$`Pr(>F)`[1])
}

cat("\nTesting", nrow(pan), "genes...\n")
tt <- t(sapply(pan$ensembl, test_one))
res <- bind_cols(pan, as.data.frame(tt)) %>%
  mutate(pct = 100 * sapply(ratio, function(x) mean(bg2 <= x))) %>%
  arrange(desc(ratio))

# BH within the post-hoc panel only, and only over NEW genes
res$q <- NA
i <- res$status == "New (post-hoc)"
res$q[i] <- p.adjust(res$p[i], method = "BH")

cat("\n===== TOP 20 BY VARIANCE RATIO =====\n")
print(res %>% dplyr::select(SYMBOL, group, status, ratio, pct, p, q) %>%
        head(20), digits = 3, row.names = FALSE)


# ---- 4. Group summaries, EXCLUDING frozen-set members ---------------

summ <- res %>%
  filter(status == "New (post-hoc)") %>%
  group_by(group) %>%
  summarise(n = n(), median_ratio = median(ratio),
            n_above_99 = sum(ratio >= quantile(bg2, .99)),
            min_q = suppressWarnings(min(q, na.rm = TRUE)),
            .groups = "drop") %>%
  arrange(desc(median_ratio))

cat("\n===== BY PATHWAY POSITION (post-hoc genes only) =====\n")
print(as.data.frame(summ), digits = 3, row.names = FALSE)

cat("\nBackground median:", round(median(bg2), 3),
    "  99th percentile:", round(quantile(bg2, .99), 3), "\n")

cat("\nHOW TO READ THIS\n")
cat("  If the fusion groups (Rab/effectors, Tethering, SNAREs, Lipid\n")
cat("  enzymes) sit above the upstream groups (signalling, initiation,\n")
cat("  conjugation), that is consistent with variance concentrating at\n")
cat("  the fusion step. It is an observation, not a test: the grouping\n")
cat("  is post-hoc and the comparison was not pre-specified.\n")
cat("\n  If nothing separates, say so. The primary result does not\n")
cat("  depend on this panel.\n")

write.csv(res, file.path(FIGS, "extended_autophagy_panel.csv"),
          row.names = FALSE)


# ---- 5. Figure ------------------------------------------------------

ord <- summ$group[order(summ$median_ratio)]
res$group <- factor(res$group,
                    levels = c(ord, setdiff(unique(res$group), ord)))

lab <- res %>% filter(ratio >= quantile(bg2, .98) |
                      SYMBOL %in% c("STX7","VPS33A","PLEKHM1","RAB7A"))

fig <- ggplot(res, aes(ratio, group)) +
  annotate("rect", xmin = quantile(bg2, .25), xmax = quantile(bg2, .75),
           ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_vline(xintercept = median(bg2), colour = "grey55", linewidth = .3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = quantile(bg2, .99), linetype = "dotted",
             colour = "grey20") +
  geom_jitter(aes(colour = status, shape = status),
              height = .18, size = 1.9, alpha = .9) +
  ggrepel::geom_text_repel(data = lab, aes(label = SYMBOL), size = 2.3,
                           seed = 1, max.overlaps = 25,
                           segment.size = .2, colour = "grey25") +
  scale_colour_manual(values = c(`New (post-hoc)` = "#5c6bc0",
                                 `In frozen crinophagy set` = "#1b5e20",
                                 `In frozen comparator set` = "#90a4ae"),
                      name = NULL) +
  scale_shape_manual(values = c(`New (post-hoc)` = 16,
                                `In frozen crinophagy set` = 17,
                                `In frozen comparator set` = 15),
                     name = NULL) +
  labs(x = "SD ratio (T2D / non-diabetic)", y = NULL,
       title = "Exploratory: variance across the autophagy-lysosome system",
       subtitle = paste0("Post-hoc panel, fully adjusted model. Grey band: ",
                         "interquartile range of the genome-wide background."),
       caption = paste(
         "POST-HOC. This panel was assembled after the primary analysis and is not part of the",
         "\npre-registered design. Triangles mark genes already in the frozen crinophagy set;",
         "\nthese are shown for position but excluded from all group summaries, since a comparator",
         "\nsharing members with the test set cannot serve as a comparator. Group medians and",
         "\nFDR values are computed over post-hoc genes only.")) +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 8),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        plot.caption = element_text(size = 6.6, colour = "grey40", hjust = 0),
        legend.position = "top")

ggsave(file.path(FIGS, "FigS_extended_autophagy.pdf"), fig,
       width = 8.5, height = 7)

cat("\nWrote FigS_extended_autophagy.pdf and extended_autophagy_panel.csv\n")
