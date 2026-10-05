# =====================================================================
# EXPLORATORY: donor-side vs acceptor-side variance
#
# STATUS: POST-HOC. Assembled after the primary analysis and after the
# GABARAP result was known. Not part of the pre-registered design. Does
# not replace the frozen 14-gene comparator. No p-value from this panel
# is presented alongside the pre-registered tests as equivalent.
#
# The frozen results stand unchanged:
#   crinophagy  31 genes, median 1.178, permutation P = 0.036
#   comparator  14 genes, median 1.036, permutation P = 0.72
#
# THE QUESTION
#
#   Crinophagy and macroautophagy converge on the same lysosome and
#   share its acceptor machinery (STX7, STX8, VTI1B, HOPS, PLEKHM1,
#   RAB7A). What differs is the donor organelle: a secretory granule
#   carrying VAMP2 with Munc13-4 and Rab27A, versus an autophagosome
#   carrying STX17/SNAP29. So: is elevated inter-donor variance a
#   property of the shared acceptor side, or of the granule donor side
#   specifically?
#
#   Munc13-4 at the 99.4th percentile and Munc13-2 at the 5.9th is a
#   hint that the answer is donor side. This asks it systematically.
#
# OVERLAP HANDLING
#
#   ~27 genes are already in the frozen crinophagy set. They are
#   flagged, plotted for position, and EXCLUDED from every tier
#   summary and every FDR calculation, because a comparison that
#   includes the test genes is not a comparison.
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

INCLUDE_TIER4 <- TRUE     # upstream signalling: regulatory, not mechanical


# ---- 1. Panel definition --------------------------------------------

tiers <- list(

  # ---------- TIER 1: autophagosome donor side ----------
  `Autophagosome: initiation` =
    c("ULK1","ULK2","ATG13","ATG101","RB1CC1"),
  `Autophagosome: nucleation` =
    c("PIK3C3","BECN1","ATG14","UVRAG","NRBF2","RUBCN","RUBCNL"),
  `Autophagosome: early membrane` =
    c("ATG9A","ATG2A","ATG2B","TMEM41B","VMP1","WIPI2","ZFYVE1",
      "VPS13A","VPS13C"),
  `Autophagosome: conjugation` =
    c("ATG3","ATG5","ATG7","ATG10","ATG12","ATG16L1","ATG4B"),
  `Autophagosome: ATG8 family` =
    c("MAP1LC3A","MAP1LC3B","GABARAP","GABARAPL1","GABARAPL2"),
  `Autophagosome: cargo receptors` =
    c("SQSTM1","NBR1","CALCOCO2","TAX1BP1","BNIP3L","OPTN"),
  `Autophagosome: fusion SNAREs` =
    c("STX17","SNAP29","YKT6","SNAP47","VAMP8","SCFD1"),
  `Autophagosome: ESCRT / closure` =
    c("VPS4A","VPS4B","VPS37A","CHMP2A","CHMP4B","VPS25","VPS28",
      "TSG101"),

  # ---------- TIER 2: granule donor side ----------
  `Granule: Munc13 priming` =
    c("UNC13D","UNC13A","UNC13B"),
  `Granule: Rab27A effectors` =
    c("RAB27A","SYTL1","SYTL2","SYTL4","RPH3AL","MLPH","MYO5A",
      "EXPH5","MADD"),
  `Granule: Rab3 axis` =
    c("RAB3A","RAB3B","RAB3D","RAB3GAP1","RAB3GAP2","RAB3IP"),
  `Granule: SNAREs` =
    c("VAMP2","VAMP4","VAMP7","VAMP3"),
  `Granule: identity and cargo` =
    c("PCLO","CHGA","CHGB","PCSK1","PCSK2","CPE","INS"),

  # ---------- TIER 3: shared lysosomal acceptor ----------
  `Lysosome: SNAREs` =
    c("STX7","STX8","VTI1B"),
  `Lysosome: HOPS` =
    c("VPS11","VPS16","VPS18","VPS33A","VPS39","VPS41"),
  `Lysosome: tethering` =
    c("PLEKHM1","EPG5","TECPR1","BIRC6"),
  `Lysosome: Rab and effectors` =
    c("RAB7A","RAB2A","RAB21","RAB24","RAB33B","RAB39A","ARL8A",
      "ARL8B","RILP","FYCO1","OSBPL1A","TBC1D25"),
  `Lysosome: membrane and lumen` =
    c("LAMP1","LAMP2","CTSB","CTSD","CTSL","ATP6V0D1","ATP6V1A",
      "TFEB","PIKFYVE","MCOLN1","NPC1"),
  `Lysosome: lipid enzymes` =
    c("INPP5E","PI4K2A","PIP5K1C")
)

if (INCLUDE_TIER4) {
  tiers[["Upstream signalling"]] <-
    c("MTOR","RPTOR","PRKAA1","PRKAA2","TBK1","STK4","IKBKB","STING1",
      "SIRT3","SIRT4","SIRT5","SIRT6","KAT2B","CREBBP","HDAC2")
}

tier_of <- function(g) {
  if (grepl("^Autophagosome", g)) "1. Autophagosome donor"
  else if (grepl("^Granule", g)) "2. Granule donor"
  else if (grepl("^Lysosome", g)) "3. Shared lysosomal acceptor"
  else "4. Upstream signalling"
}

pan <- bind_rows(lapply(names(tiers), function(g)
  data.frame(SYMBOL = tiers[[g]], subgroup = g,
             tier = tier_of(g), stringsAsFactors = FALSE)
)) %>% distinct(SYMBOL, .keep_all = TRUE)

cat("Genes proposed:", nrow(pan), "\n")
print(table(pan$tier))


# ---- 2. Flag frozen-set overlap -------------------------------------

frozen <- read.delim("frozen_gene_sets_v2_clean.txt",
                     stringsAsFactors = FALSE)
cr <- frozen$SYMBOL[frozen$set == "crinophagy"]
cm <- frozen$SYMBOL[frozen$set == "comparator"]

pan$status <- ifelse(pan$SYMBOL %in% cr, "Frozen crinophagy set",
              ifelse(pan$SYMBOL %in% cm, "Frozen comparator set",
                     "New (post-hoc)"))

cat("\n===== FROZEN-SET OVERLAP =====\n")
print(table(pan$status))
cat("\nExcluded from all tier summaries:\n  ",
    paste(pan$SYMBOL[pan$status != "New (post-hoc)"], collapse = ", "),
    "\n")


# ---- 3. Map and test ------------------------------------------------

mp <- AnnotationDbi::select(org.Hs.eg.db, keys = pan$SYMBOL,
                            keytype = "SYMBOL", columns = "ENSEMBL")
mp <- mp[!is.na(mp$ENSEMBL) & mp$ENSEMBL %in% rownames(mat), ]
mp <- mp[!duplicated(mp$SYMBOL), ]

absent <- setdiff(pan$SYMBOL, mp$SYMBOL)
pan <- pan %>% inner_join(mp, by = "SYMBOL") %>% rename(ensembl = ENSEMBL)

cat("\nMeasured:", nrow(pan), " Not measured:", length(absent), "\n")
if (length(absent)) cat("  ", paste(absent, collapse = ", "), "\n")

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

cat("\nTesting", nrow(pan), "genes (bootstrap 2000 each)...\n")
tt <- t(sapply(pan$ensembl, test_one))
res <- bind_cols(pan, as.data.frame(tt)) %>%
  mutate(pct = 100 * sapply(ratio, function(x) mean(bg2 <= x)))

res$q <- NA
i <- res$status == "New (post-hoc)"
res$q[i] <- p.adjust(res$p[i], method = "BH")

res <- res %>% arrange(desc(ratio))

cat("\n===== TOP 25 =====\n")
print(res %>% dplyr::select(SYMBOL, subgroup, status, ratio, pct, p, q) %>%
        head(25), digits = 3, row.names = FALSE)


# ---- 4. Tier and subgroup summaries (post-hoc genes only) -----------

p99 <- quantile(bg2, .99)

tier_summ <- res %>% filter(status == "New (post-hoc)") %>%
  group_by(tier) %>%
  summarise(n = n(), median = median(ratio),
            above99 = sum(ratio >= p99), .groups = "drop") %>%
  arrange(desc(median))

sub_summ <- res %>% filter(status == "New (post-hoc)") %>%
  group_by(tier, subgroup) %>%
  summarise(n = n(), median = median(ratio),
            above99 = sum(ratio >= p99), .groups = "drop") %>%
  arrange(desc(median))

cat("\n===== BY TIER (post-hoc genes only) =====\n")
print(as.data.frame(tier_summ), digits = 3, row.names = FALSE)

cat("\n===== BY SUBGROUP =====\n")
print(as.data.frame(sub_summ), digits = 3, row.names = FALSE)

cat("\nBackground median:", round(median(bg2), 3),
    " 99th pct:", round(p99, 3), "\n")

cat("\nHOW TO READ THIS\n")
cat("  Donor-side hypothesis: granule tier above the shared lysosomal\n")
cat("  tier. Note that the granule tier loses most of its members to\n")
cat("  the frozen-set exclusion, so its post-hoc n is small; the\n")
cat("  comparison it supports is weak and descriptive.\n")
cat("  If the tiers do not separate, report that. The primary result\n")
cat("  does not depend on this panel.\n")

write.csv(res, file.path(FIGS, "donor_acceptor_panel.csv"),
          row.names = FALSE)


# ---- 5. Figure ------------------------------------------------------

ordsub <- sub_summ$subgroup[order(sub_summ$median)]
extra <- setdiff(unique(res$subgroup), ordsub)
res$subgroup <- factor(res$subgroup, levels = c(ordsub, extra))

lab <- res %>%
  filter(ratio >= quantile(bg2, .975) | ratio <= quantile(bg2, .03) |
         SYMBOL %in% c("UNC13D","STX7","RAB27A","VAMP2","RAB7A"))

fig <- ggplot(res, aes(ratio, subgroup)) +
  annotate("rect", xmin = quantile(bg2, .25), xmax = quantile(bg2, .75),
           ymin = -Inf, ymax = Inf, fill = "grey93") +
  geom_vline(xintercept = median(bg2), colour = "grey55", linewidth = .3) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") +
  geom_vline(xintercept = p99, linetype = "dotted", colour = "grey20") +
  geom_jitter(aes(colour = status, shape = status),
              height = .16, size = 1.8, alpha = .9) +
  ggrepel::geom_text_repel(data = lab, aes(label = SYMBOL), size = 2.2,
                           seed = 1, max.overlaps = 30,
                           segment.size = .2, colour = "grey25") +
  facet_grid(tier ~ ., scales = "free_y", space = "free_y",
             switch = "y") +
  scale_colour_manual(values = c(`New (post-hoc)` = "#5c6bc0",
                                 `Frozen crinophagy set` = "#1b5e20",
                                 `Frozen comparator set` = "#90a4ae"),
                      name = NULL) +
  scale_shape_manual(values = c(`New (post-hoc)` = 16,
                                `Frozen crinophagy set` = 17,
                                `Frozen comparator set` = 15),
                     name = NULL) +
  labs(x = "SD ratio (T2D / non-diabetic)", y = NULL,
       title = "Exploratory: donor-side and acceptor-side machinery",
       subtitle = "Fully adjusted model. Grey band: interquartile range of the genome-wide background.",
       caption = paste(
         "POST-HOC PANEL, assembled after the primary analysis; not pre-registered. Crinophagy and",
         "\nmacroautophagy share the lysosomal acceptor machinery (tier 3) but differ in the donor",
         "\norganelle (tiers 1 and 2). Triangles mark frozen crinophagy-set members: shown for",
         "\nposition, excluded from all tier medians and FDR values, since a comparison containing",
         "\nthe test genes is not a comparison.")) +
  theme_minimal(base_size = 9) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 7.5),
        strip.text.y.left = element_text(angle = 0, size = 8,
                                         face = "bold", colour = "grey25"),
        strip.placement = "outside",
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8, colour = "grey35"),
        plot.caption = element_text(size = 6.4, colour = "grey40", hjust = 0),
        legend.position = "top")

ggsave(file.path(FIGS, "FigS_donor_acceptor.pdf"), fig,
       width = 8.5, height = 10)

cat("\nWrote FigS_donor_acceptor.pdf and donor_acceptor_panel.csv\n")
