# 71_harmonized_raw_31_crinophagy_audit.R
# IJMS revision: investigate five fully measured raw features absent from the
# processed HumanIslets matrix; reanalyse the entire frozen 31-gene crinophagy
# panel using ONE explicitly defined raw-log2/median-centred scale.
#
# READ ONLY with respect to originals: writes 41-44 under reviewer_validation_output.
# WARNING: this is a NEW exploratory sensitivity analysis, NOT a reconstruction
# of the provider's exact filtering, imputation, technical-replicate handling,
# nor a validated covariate-adjusted dispersion hypothesis test.

ROOT <- "."
OUT <- file.path(ROOT, "reviewer_validation_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

needed_files <- c(
  file.path(ROOT, "proteomics", "unproc_prot.rds"),
  file.path(ROOT, "proteomics", "proc_prot.rds"),
  file.path(ROOT, "donor.csv"),
  file.path(OUT, "15_raw_uniprot_entrez_mapping.csv"),
  file.path(OUT, "16_raw_and_processed_panel_coverage.csv"),
  file.path(ROOT, "per_protein_variance.csv")
)
if (any(!file.exists(needed_files)))
  stop("Required inputs missing:\n", paste(needed_files[!file.exists(needed_files)], collapse="\n"))

read_wide <- function(path, idcol) {
  a <- readRDS(path)
  if (!is.data.frame(a) || !(idcol %in% colnames(a)))
    stop("Expected data.frame with ", idcol, ": ", path)
  ids <- as.character(a[[idcol]])
  if (anyNA(ids) || anyDuplicated(ids) || any(!nzchar(ids))) stop("Invalid feature IDs: ", path)
  x <- as.matrix(a[,setdiff(colnames(a),idcol),drop=FALSE])
  storage.mode(x) <- "numeric"
  rownames(x) <- ids
  if (is.null(colnames(x)) || anyDuplicated(colnames(x))) stop("Invalid donor IDs in ", path)
  x
}

RAW <- read_wide(needed_files[1], "feature_id")
PRO <- read_wide(needed_files[2], "gene_id")
mp <- read.csv(needed_files[4], colClasses="character", check.names=FALSE)
cv <- read.csv(needed_files[5], stringsAsFactors=FALSE, check.names=FALSE)
cv <- cv[cv$set=="crinophagy",,drop=FALSE]
if (nrow(cv)!=31L || anyDuplicated(cv$gene))
  stop("Expected 31 unique crinophagy genes in file 16")
if (!all(c("feature_id","mapped_entrez","status") %in% names(mp)))
  stop("File 15 does not contain required mapping columns")
if (!all(c("gene","entrez","annotation_status","processed_rows") %in% names(cv)))
  stop("File 16 lacks required coverage columns")

D <- read.csv(needed_files[3], stringsAsFactors=FALSE, check.names=FALSE)
COV <- c("donorage", "donorsex", "bodymassindex",
         "predistributionculturetime", "percentieqrecoverypostculture")
req <- c("record_id", "diagnosis_computed", COV)
if (!all(req %in% names(D))) stop("Missing metadata: ", paste(setdiff(req,names(D)),collapse=", "))
if (anyDuplicated(D$record_id)) stop("Duplicate donor metadata identifiers")
D$group <- unname(c("None"="ND", "Pre.T2D"="preT2D", "Type2"="T2D") [
  trimws(as.character(D$diagnosis_computed))])
D <- D[complete.cases(D[,COV,drop=FALSE]) & !is.na(D$group) &
         D$record_id %in% intersect(colnames(RAW),colnames(PRO)), ,drop=FALSE]
D$group <- factor(D$group,levels=c("ND","preT2D","T2D"))
D$donorsex <- factor(D$donorsex)
if (!identical(unname(as.integer(table(D$group))),c(87L,20L,16L)))
  stop("Expected 87 ND, 20 preT2D, 16 T2D covariate-complete donors")
ids <- as.character(D$record_id)

# Match Script 67's exploratory normalization for EVERY raw feature:
# log2(raw intensity) then subtract EACH DONOR'S feature median on that scale
# and add the across-donor grand median. No imputation.
X <- RAW
X[!is.finite(X) | X<=0] <- NA_real_
L <- log2(X)
med <- apply(L,2,median,na.rm=TRUE)
if (any(!is.finite(med))) stop("A donor lacks finite positive raw signals")
LN <- sweep(L,2,med,"-") + median(med)

# Raw feature variance rank over all 129 donors, not merely the 123 fitted.
# This is an approximation of possible low-variance filtering, NOT a
# verified replication of the provider's processed workflow.
nraw <- rowSums(is.finite(LN))
eligible <- nraw >= ceiling(0.5*ncol(LN))
feature_var <- apply(LN,1,var,na.rm=TRUE)
reference <- feature_var[eligible & is.finite(feature_var)]
if (length(reference)<100) stop("Too few reference features for variance ranking")
var_bottom_20 <- as.numeric(quantile(reference,0.20,na.rm=TRUE,names=FALSE))

fit_one <- function(y) {
  ok <- is.finite(y)
  dd <- D[ok,,drop=FALSE]
  dd$y <- as.numeric(y[ok])
  nd <- dd$group=="ND"; td <- dd$group=="T2D"
  n_nd <- sum(nd); n_td <- sum(td)
  out <- c(n_model=nrow(dd),n_ND=n_nd,n_preT2D=sum(dd$group=="preT2D"),
           n_T2D=n_td, sd_ratio=NA_real_, bf_p_approx=NA_real_,
           sd_ratio_group_mean=NA_real_)
  if (n_nd<5L || n_td<5L || sum(ok)<25L) return(out)
  form <- y ~ donorage + donorsex + bodymassindex +
    predistributionculturetime + percentieqrecoverypostculture
  design <- model.matrix(form,data=dd)
  if (qr(design)$rank != ncol(design)) return(out)
  r <- residuals(lm(form,data=dd))
  sn <- sd(r[nd]); st <- sd(r[td])
  if (!is.finite(sn) || sn<=0 || !is.finite(st)) return(out)
  out["sd_ratio"] <- st/sn
  zz <- c(abs(r[nd]-median(r[nd])),abs(r[td]-median(r[td])))
  lab <- factor(c(rep("ND",n_nd),rep("T2D",n_td)),levels=c("ND","T2D"))
  out["bf_p_approx"] <- tryCatch(
    as.numeric(anova(lm(zz~lab))[["Pr(>F)"]][1]),error=function(e) NA_real_)
  group_form <- update(form,. ~ . + group)
  gX <- model.matrix(group_form,data=dd)
  if (qr(gX)$rank == ncol(gX)) {
    rg <- residuals(lm(group_form,data=dd))
    sg <- sd(rg[nd]);
    if (is.finite(sg) && sg>0) out["sd_ratio_group_mean"] <- sd(rg[td])/sg
  }
  out
}

accession_for <- function(ent) {
  if (is.na(ent) || !nzchar(ent)) return(character())
  idx <- !is.na(mp$mapped_entrez) & mp$status=="UNIQUE" & mp$mapped_entrez==ent
  intersect(unique(mp$feature_id[idx]),rownames(LN))
}

per_gene <- vector("list",nrow(cv))
per_feature <- list()
for (i in seq_len(nrow(cv))) {
  g <- as.character(cv$gene[i]); ent <- as.character(cv$entrez[i])
  if (is.na(cv$entrez[i]) || cv$annotation_status[i]!="unique_ensembl_entrez")
    ent <- NA_character_
  acc <- accession_for(ent)
  y <- rep(NA_real_,nrow(D))
  raw_rule <- "none_mapping_unresolved"
  if (length(acc)==1L) {
    y <- as.numeric(LN[acc,ids]); raw_rule <- "single_uniquely_mapped_accession"
  } else if (length(acc)>1L) {
    # Explicit exploratory aggregation; not asserted to match provider.
    tmp <- LN[acc,ids,drop=FALSE]
    y <- colMeans(tmp,na.rm=TRUE); y[!is.finite(y)] <- NA_real_
    raw_rule <- "multiple_accessions_mean_log2_EXPLORATORY"
  }
  z <- fit_one(y)
  per_gene[[i]] <- data.frame(gene=g,entrez=ent,
    n_raw_accessions=length(acc),raw_accessions=paste(acc,collapse=";"),
    raw_aggregation=raw_rule,processed_rows=cv$processed_rows[i],
    n_model=as.integer(z["n_model"]),n_ND=as.integer(z["n_ND"]),
    n_preT2D=as.integer(z["n_preT2D"]),n_T2D=as.integer(z["n_T2D"]),
    sd_ratio=as.numeric(z["sd_ratio"]),bf_p_approx=as.numeric(z["bf_p_approx"]),
    sd_ratio_group_mean=as.numeric(z["sd_ratio_group_mean"]),
    stringsAsFactors=FALSE)
  if (length(acc)) for (a in acc) {
    v <- feature_var[a]
    per_feature[[length(per_feature)+1L]] <- data.frame(
      gene=g,entrez=ent,accession=a,
      n_positive_finite_all_raw_donors=as.integer(nraw[a]),
      positive_finite_fraction_all_raw_donors=as.numeric(nraw[a]/ncol(LN)),
      log2_median_centered_variance=as.numeric(v),
      variance_percentile_from_low=as.numeric(mean(reference<=v)),
      at_or_below_20pct_reference=as.logical(v<=var_bottom_20),
      processed_rows=cv$processed_rows[i],stringsAsFactors=FALSE)
  }
}
result <- do.call(rbind,per_gene)
# Each member of frozen family receives an explicit row. Unevaluable cases
# receive p=NA; for the 31-test conservative sensitivity, substitute p=1.
result$q_BH_evaluable <- NA_real_
idx <- which(is.finite(result$bf_p_approx))
if (length(idx)) result$q_BH_evaluable[idx] <- p.adjust(result$bf_p_approx[idx],method="BH")
p_for31 <- result$bf_p_approx
p_for31[!is.finite(p_for31)] <- 1
result$q_BH_31_p1_for_unevaluable <- p.adjust(p_for31,method="BH")
result$discovery_orig_tested_processed <- result$gene %in% c(
  read.csv(needed_files[6],stringsAsFactors=FALSE)$symbol)
write.csv(result,file.path(OUT,"41_raw31_variance_and_BH_sensitivity.csv"),row.names=FALSE)

feature_result <- if(length(per_feature)) do.call(rbind,per_feature) else
  data.frame(gene=character(),entrez=character(),accession=character())
write.csv(feature_result,file.path(OUT,"42_raw31_feature_variance_rank.csv"),row.names=FALSE)

five <- c("STX7","VPS33A","VPS11","VAMP7","VPS16")
cat("\n===== ALL 31 RAW-PANEL VARIANCE TESTS (EXPLORATORY) =====\n")
print(result[,c("gene","n_raw_accessions","n_ND","n_T2D","sd_ratio",
               "bf_p_approx","q_BH_evaluable","sd_ratio_group_mean")],row.names=FALSE,digits=4)
cat("\n===== FIVE EXCLUDED PROCESSED FEATURES =====\n")
print(result[result$gene %in% five,c("gene","raw_accessions","n_ND","n_T2D",
                                      "sd_ratio","bf_p_approx","q_BH_evaluable")],
      row.names=FALSE,digits=5)
cat("\n===== FIVE FEATURE VARIANCE RANKS (LOW TO HIGH) =====\n")
print(feature_result[feature_result$gene %in% five,
  c("gene","accession","n_positive_finite_all_raw_donors",
    "variance_percentile_from_low","at_or_below_20pct_reference")],row.names=FALSE,digits=4)
cat(sprintf("\nFeatures eligible for raw variance reference: %d. Approx bottom-quintile variance threshold: %.6f\n",
            length(reference),var_bottom_20))
cat(sprintf("Raw 31-panel evaluable gene tests: %d of 31.\n",length(idx)))

notes <- c(
  "This analysis used the ORIGINAL raw UniProt matrix, not the processed matrix.",
  "Steps: positive finite raw values -> log2 -> subtract each donor's across-feature median -> add overall median; no imputation.",
  "Multi-accession genes were averaged on transformed log2 scale: exploratory, not provider-verified aggregation.",
  "Covariate-only nuisance model fitted over available ND/preT2D/T2D donors; SD ratios computed ND vs T2D.",
  "Group-aware nuisance model is a separate descriptive sensitivity.",
  "Brown-Forsythe P values from median-centred ABS residuals are approximate, not validated conditional covariate inference.",
  "BH within evaluable raw crinophagy panel is exploratory; q_BH_31_p1_for_unevaluable imputes p=1 as a separate conservative family-size sensitivity.",
  "Raw variance ranks across >=50%-observed features are a proxy for potential low-variance filtering, not proof of provider exclusion.",
  "Provider filtering may have occurred after technical replicate handling, median normalization, missForest imputation, and identifier aggregation.",
  "A zero raw NA fraction does NOT mean a gene necessarily passes all processing filters; the five raw-complete genes could be low-variance or removed at other stages.",
  "Do not replace the historical 26-processed-protein q values by these new raw sensitivity q values without explicitly stating the new analysis.",
  "Do not treat q<0.05 here as confirmatory if source version, batch, leverage, missingness, and resampling calibration remain unresolved."
)
writeLines(notes,file.path(OUT,"43_raw31_interpretation_cautions.txt"))
write.csv(data.frame(
  item=c("raw_features","processed_features","donors_raw","donors_complete_cov",
         "eligible_feature_variance_reference","raw_feature_bottom20_var_threshold",
         "crinophagy_genes","raw_genes_evaluable"),
  value=c(nrow(RAW),nrow(PRO),ncol(RAW),nrow(D),length(reference),var_bottom_20,31L,length(idx))),
  file.path(OUT,"44_raw31_analysis_manifest.csv"),row.names=FALSE)
cat("\nCreated files 41, 42, 43, 44 under reviewer_validation_output.\n")
