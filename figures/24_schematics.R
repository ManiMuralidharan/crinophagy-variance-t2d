# =====================================================================
# Two schematics as editable vector PDFs
#
#   Fig 1B  Crinophagy machinery, with frozen-set members marked
#   Fig 3D  Munc13 paralogue divergence: two destinations, one family
#
# EDITABILITY. Output is written with grDevices::pdf(), so every
# element is a vector object with live text. Opening either file in
# Illustrator or Inkscape gives editable paths and editable type; no
# rasterisation, no outlined fonts. If your Illustrator shows text as
# outlines, set useDingbats = FALSE (done below) and check that the
# font is embedded rather than substituted.
#
# WHY DRAWN RATHER THAN DIAGRAMMED. These use base grid graphics so
# the geometry is explicit and adjustable. Nothing here is data; both
# are illustrations of published mechanism (ours) and of the analysis
# design.
# =====================================================================

setwd("/Users/mmani3/Downloads/Human Islets")
FIGS <- "figures"
dir.create(FIGS, showWarnings = FALSE)

library(grid)

# palette, matched to the data figures
GREEN  <- "#1b5e20"   # crinophagy / frozen set
GREY   <- "#90a4ae"   # comparator / not in set
ORANGE <- "#ef6c00"   # highlighted
LIGHT  <- "#e8f5e9"
SAND   <- "#fff8e1"
INK    <- "#212121"


# =====================================================================
# FIGURE 1B - crinophagy machinery
# =====================================================================

pdf(file.path(FIGS, "Fig1B_crinophagy_schematic.pdf"),
    width = 7, height = 4.6, useDingbats = FALSE)

grid.newpage()
pushViewport(viewport(width = .94, height = .94))

grid.text("Crinophagy: secretory granule to lysosome",
          x = .02, y = .96, just = "left",
          gp = gpar(fontsize = 11, fontface = "bold", col = INK))
grid.text("Bold: member of the pre-registered 31-protein set",
          x = .02, y = .905, just = "left",
          gp = gpar(fontsize = 7.5, col = "grey40"))

# --- lysosome ---
grid.circle(x = .70, y = .50, r = .20,
            gp = gpar(fill = SAND, col = "grey55", lwd = 1.1))
grid.text("Lysosome", x = .70, y = .50,
          gp = gpar(fontsize = 10, col = "grey30"))

# --- granule, docked ---
grid.circle(x = .38, y = .50, r = .095,
            gp = gpar(fill = LIGHT, col = GREEN, lwd = 1.3))
grid.text("Insulin\ngranule", x = .38, y = .50,
          gp = gpar(fontsize = 8, col = GREEN, lineheight = .95))

# --- granule approaching ---
grid.circle(x = .14, y = .50, r = .075,
            gp = gpar(fill = LIGHT, col = "grey60", lwd = 1))
grid.lines(x = c(.225, .275), y = c(.50, .50),
           arrow = arrow(length = unit(.06, "inches"), type = "closed"),
           gp = gpar(col = "grey50", lwd = 1.1, fill = "grey50"))

# --- granule-side machinery ---
gran <- c("RAB27A", "Munc13-4", "VAMP2")
for (i in seq_along(gran)) {
  y <- .74 - (i - 1) * .075
  grid.text(gran[i], x = .30, y = y, just = "right",
            gp = gpar(fontsize = 8.5, fontface = "bold", col = GREEN))
}
grid.lines(x = c(.315, .355), y = c(.665, .565),
           gp = gpar(col = GREEN, lwd = .7, lty = "dotted"))
grid.text("granule side", x = .30, y = .795, just = "right",
          gp = gpar(fontsize = 7, col = "grey45", fontface = "italic"))

# --- lysosome-side machinery ---
lys <- c("STX7", "STX8", "VTI1B", "HOPS (VPS11/16/33A)", "PLEKHM1")
for (i in seq_along(lys)) {
  y <- .78 - (i - 1) * .068
  bold <- lys[i] %in% c("STX7", "STX8", "HOPS (VPS11/16/33A)")
  grid.text(lys[i], x = .885, y = y, just = "left",
            gp = gpar(fontsize = 8.5,
                      fontface = if (bold) "bold" else "plain",
                      col = if (bold) GREEN else GREY))
}
grid.lines(x = c(.875, .805), y = c(.70, .60),
           gp = gpar(col = GREEN, lwd = .7, lty = "dotted"))
grid.text("lysosome side", x = .885, y = .845, just = "left",
          gp = gpar(fontsize = 7, col = "grey45", fontface = "italic"))

# --- fusion interface ---
grid.lines(x = c(.475, .50), y = c(.50, .50),
           gp = gpar(col = ORANGE, lwd = 2))
grid.text("Ca2+", x = .4875, y = .565,
          gp = gpar(fontsize = 8, col = ORANGE, fontface = "bold"))

# --- degradation output ---
grid.text("cargo degraded", x = .70, y = .295,
          gp = gpar(fontsize = 8, col = "grey35"))
grid.text("CTSB  CTSD  CTSL", x = .70, y = .245,
          gp = gpar(fontsize = 8.5, fontface = "bold", col = GREEN))
grid.text("membrane recycled", x = .70, y = .195,
          gp = gpar(fontsize = 8, col = "grey35"))

# --- footnote ---
grid.text("Machinery as defined in our previous work (ref).",
          x = .02, y = .045, just = "left",
          gp = gpar(fontsize = 7, col = "grey45"))

popViewport()
dev.off()

cat("Wrote Fig1B_crinophagy_schematic.pdf\n")


# =====================================================================
# FIGURE 3D - Munc13 paralogue divergence
# =====================================================================

pdf(file.path(FIGS, "Fig3D_munc13_schematic.pdf"),
    width = 6.6, height = 4.2, useDingbats = FALSE)

grid.newpage()
pushViewport(viewport(width = .94, height = .94))

grid.text("One granule, two destinations, three paralogues",
          x = .02, y = .955, just = "left",
          gp = gpar(fontsize = 11, fontface = "bold", col = INK))

# --- central granule ---
grid.circle(x = .20, y = .50, r = .085,
            gp = gpar(fill = LIGHT, col = "grey45", lwd = 1.2))
grid.text("Insulin\ngranule", x = .20, y = .50,
          gp = gpar(fontsize = 8, col = "grey25", lineheight = .95))

# --- upper branch: exocytosis ---
grid.curve(x1 = .29, y1 = .545, x2 = .55, y2 = .755,
           curvature = -.3, square = FALSE,
           arrow = arrow(length = unit(.07, "inches"), type = "closed"),
           gp = gpar(col = GREY, lwd = 1.5, fill = GREY))
grid.rect(x = .755, y = .755, width = .30, height = .085,
          gp = gpar(fill = "grey95", col = GREY, lwd = 1))
grid.text("Plasma membrane", x = .755, y = .775,
          gp = gpar(fontsize = 9, col = "grey25"))
grid.text("exocytosis", x = .755, y = .733,
          gp = gpar(fontsize = 7.5, col = "grey45", fontface = "italic"))

grid.text("Munc13-1  Munc13-2", x = .42, y = .675,
          gp = gpar(fontsize = 8.5, col = GREY, fontface = "bold"))
grid.text("51.5th and 5.9th percentile", x = .42, y = .632,
          gp = gpar(fontsize = 7.5, col = "grey50"))

# --- lower branch: crinophagy ---
grid.curve(x1 = .29, y1 = .455, x2 = .55, y2 = .245,
           curvature = .3, square = FALSE,
           arrow = arrow(length = unit(.07, "inches"), type = "closed"),
           gp = gpar(col = GREEN, lwd = 1.8, fill = GREEN))
grid.circle(x = .755, y = .245, r = .105,
            gp = gpar(fill = SAND, col = GREEN, lwd = 1.3))
grid.text("Lysosome", x = .755, y = .265,
          gp = gpar(fontsize = 9, col = INK))
grid.text("crinophagy", x = .755, y = .222,
          gp = gpar(fontsize = 7.5, col = GREEN, fontface = "italic"))

grid.text("Munc13-4", x = .42, y = .325,
          gp = gpar(fontsize = 9.5, col = GREEN, fontface = "bold"))
grid.text("99.4th percentile", x = .42, y = .282,
          gp = gpar(fontsize = 7.5, col = GREEN))

# --- caption ---
grid.text(paste("Percentiles are inter-donor variance ratios (T2D vs non-diabetic)",
                "against a 20,255-gene background."),
          x = .02, y = .075, just = "left",
          gp = gpar(fontsize = 7, col = "grey45"))
grid.text("The paralogue directing granules to lysosomes is the one that diverges.",
          x = .02, y = .035, just = "left",
          gp = gpar(fontsize = 7, col = "grey45", fontface = "italic"))

popViewport()
dev.off()

cat("Wrote Fig3D_munc13_schematic.pdf\n")
cat("\nBoth are vector PDFs with live text. Open in Illustrator or\n")
cat("Inkscape to adjust geometry, wording or colour.\n")
