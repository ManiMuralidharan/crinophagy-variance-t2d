library(car)
set.seed(1)
library(org.Hs.eg.db); library(AnnotationDbi)
mp <- AnnotationDbi::select(org.Hs.eg.db,
        keys = c("CADPS","CADPS2","BAIAP3"),
        keytype = "SYMBOL", columns = "ENSEMBL")
mp <- mp[!is.na(mp$ENSEMBL) & mp$ENSEMBL %in% rownames(mat), ]
mp <- mp[!duplicated(mp$SYMBOL), ]
print(mp)
detach("package:org.Hs.eg.db", unload = TRUE)
detach("package:AnnotationDbi", unload = TRUE)
for (i in seq_len(nrow(mp))) {
  e <- mp$ENSEMBL[i]
  d <- pd2
  d$y <- mat[e, rownames(d)]
  d$r <- residuals(lm(y ~ batch + imm + beta, data = d))
  a <- d$r[d$group == "T2D"]; b <- d$r[d$group == "ND"]
  bt <- replicate(5000,
    sd(sample(a, length(a), TRUE)) / sd(sample(b, length(b), TRUE)))
  cat(sprintf("%-8s ratio=%.3f  CI[%.3f, %.3f]  p=%.4f  pct=%.1f\n",
      mp$SYMBOL[i], sd(a)/sd(b),
      quantile(bt,.025), quantile(bt,.975),
      leveneTest(r ~ group, data = d)$`Pr(>F)`[1],
      100*mean(bg2 <= sd(a)/sd(b))))
}
cat("\nNot measured:",
    paste(setdiff(c("CADPS","CADPS2","BAIAP3"), mp$SYMBOL), collapse=", "), "\n")
