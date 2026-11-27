# 60_reviewer_replication_validation.R
# IJMS revision: R1-M1/M2/M3; R2-M1/M2/M3/M6/M8 and related minor points.
# Tested design: replication GSE164416, 18 ND + 36 T2D (check, never assume).
# PURPOSE: reproduce multiplicity and assess covariate/model/donor sensitivity.
# None of the newly computed tests are automatically confirmatory.
# Run from any location: Rscript 60_reviewer_replication_validation.R

ROOT <- "."  # CHANGE ONLY IF NEEDED
SESSION <- file.path(ROOT, "session_full.RData")
S1FILE <- file.path(ROOT, "ijms-4686697-supplementary", "Supplementary_Tables.xlsx")
OUTPUT <- file.path(ROOT, "reviewer_validation_output")
NBOOT <- as.integer(Sys.getenv("REVIEWER_NBOOT", "1200"))
set.seed(20261113)

if (!requireNamespace("readxl", quietly = TRUE))
  stop("Package readxl is required. Run install.packages('readxl').")
stopifnot(file.exists(SESSION), file.exists(S1FILE), NBOOT >= 100L)
dir.create(OUTPUT, recursive = TRUE, showWarnings = FALSE)

save_csv <- function(x, filename) {
  write.csv(x, file.path(OUTPUT, filename), row.names = FALSE, na = "")
}
check_columns <- function(dat, cols, label) {
  missing <- setdiff(cols, colnames(dat))
  if (length(missing)) stop(label, " missing columns: ", paste(missing, collapse = ", "))
}
num <- function(x) suppressWarnings(as.numeric(x))
finite_sd_ratio <- function(a, b) {
  if (length(a) < 3 || length(b) < 3 || anyNA(c(a,b))) return(NA_real_)
  sa <- sd(a); sb <- sd(b)
  if (!is.finite(sa) || !is.finite(sb) || sb <= 0) return(NA_real_)
  sa / sb
}

cat("\n=== A. Replication panel, exact nominal and adjusted P values ===\n")
S1_all <- as.data.frame(readxl::read_excel(S1FILE, sheet="S1_pre-registered_sets", skip=2))
check_columns(S1_all, c("gene", "set", "ensembl", "p_M1", "p_M2", "sd_ratio_M2", "ci_lo_M2", "ci_hi_M2"), "S1")
# The worksheet includes STX1A as a 46th, explicitly post-hoc row.
# Preserve that provenance but exclude it from the frozen 31 + 14 panel tests.
posthoc_rows <- S1_all[!is.na(S1_all$set) &
                         !(S1_all$set %in% c("crinophagy", "comparator")),
                         c("gene", "set"), drop=FALSE]
cat("Workbook rows:", nrow(S1_all),
    "; excluded non-frozen rows:", nrow(posthoc_rows), "\n")
if (nrow(posthoc_rows)) print(posthoc_rows, row.names=FALSE)
S1 <- S1_all[!is.na(S1_all$set) &
             S1_all$set %in% c("crinophagy", "comparator"), , drop=FALSE]

S1$gene <- as.character(S1$gene)
S1$set <- as.character(S1$set)
S1$ensembl <- as.character(S1$ensembl)
stopifnot(nrow(S1) == 45L, sum(S1$set == "crinophagy") == 31L,
          sum(S1$set == "comparator") == 14L,
          !anyDuplicated(S1$gene), !anyDuplicated(S1$ensembl))
FAMILIES <- c("crinophagy", "comparator")
MULT <- do.call(rbind, lapply(c("M1", "M2"), function(model) {
  p <- num(S1[[paste0("p_", model)]])
  d <- data.frame(gene=S1$gene, ensembl=S1$ensembl, set=S1$set,
                  model=model, sd_ratio=num(S1[[paste0("sd_ratio_", model)]]),
                  ci_lo=num(S1[[paste0("ci_lo_", model)]]),
                  ci_hi=num(S1[[paste0("ci_hi_", model)]]),
                  p_nominal=p, q_joint45=p.adjust(p, method="BH"))
  d$q_own_set <- NA_real_
  for (g in FAMILIES) {
    ii <- which(d$set == g)
    d$q_own_set[ii] <- p.adjust(d$p_nominal[ii], method="BH")
  }
  # Both selected discovery candidates: one family of two (if selected before replication).
  two <- which(d$gene %in% c("UNC13D", "RAB3A"))
  if (length(two) != 2L) stop("Could not find both nominated discovery candidates")
  d$p_2candidate_holm <- NA_real_
  d$p_2candidate_bonf <- NA_real_
  d$p_2candidate_holm[two] <- p.adjust(d$p_nominal[two], method="holm")
  d$p_2candidate_bonf[two] <- p.adjust(d$p_nominal[two], method="bonferroni")
  d
}))
save_csv(MULT, "01_replication_multiplicity.csv")
print(MULT[MULT$model == "M2" & MULT$gene %in% c("UNC13D","RAB3A","STX7","GABARAP"),
           c("gene","sd_ratio","p_nominal","p_2candidate_holm","q_own_set","q_joint45")], row.names=FALSE)
cat("NOTE: Candidate prespecification must be proven by dated records, NOT by this script.\n")

cat("\n=== B. Load RNA matrix and assess batch and composition ===\n")
env <- new.env(parent=globalenv())
load(SESSION, envir=env)
for (nm in c("mat", "pd2")) if (!exists(nm, envir=env, inherits=FALSE)) stop("Missing ", nm, " in session_full.RData")
mat <- as.matrix(env$mat)
pd2 <- as.data.frame(env$pd2)
check_columns(pd2, c("group","batch","imm","beta"), "pd2")
if (is.null(rownames(pd2)) || is.null(colnames(mat))) stop("Donor IDs must be present")
if (!all(rownames(pd2) %in% colnames(mat))) stop("Donor identifiers not matched to expression matrix")
if (!all(S1$ensembl %in% rownames(mat))) stop("Panel gene Ensembl IDs absent from mat: ",
                                              paste(S1$gene[!S1$ensembl %in% rownames(mat)], collapse=", "))
if (anyNA(pd2[,c("group","batch","imm","beta")])) stop("Incomplete covariates. Resolve explicitly, do not drop silently.")
pd2$group <- factor(as.character(pd2$group), levels=c("ND","T2D"))
if (anyNA(pd2$group)) stop("Unexpected diabetes categories; check cohort selection")
pd2$batch <- factor(as.character(pd2$batch))
if (!is.numeric(pd2$imm) || !is.numeric(pd2$beta)) stop("imm and beta must be numeric")
SAMPLES <- rownames(pd2)
Y <- t(mat[S1$ensembl, SAMPLES, drop=FALSE])
colnames(Y) <- S1$gene
if (anyNA(Y) || any(!is.finite(Y))) stop("Missing/nonfinite panel values. Explicit missing-data policy needed.")
TAB <- as.data.frame.matrix(table(batch=pd2$batch, group=pd2$group))
TAB$batch <- rownames(TAB)
TAB <- TAB[, c("batch", "ND", "T2D")]
save_csv(TAB, "02_batch_status_counts.csv")
print(TAB, row.names=FALSE)
cat("All four batches should contain both classes; within-batch label exchangeability also depends on covariate structure.\n")
BAL <- do.call(rbind, lapply(c("imm", "beta"), function(v) {
  do.call(rbind, lapply(levels(pd2$batch), function(b) {
    do.call(rbind, lapply(levels(pd2$group), function(g) {
      x <- pd2[[v]][pd2$batch == b & pd2$group == g]
      data.frame(variable=v, batch=b, group=g, n=length(x), mean=mean(x), sd=sd(x))
    }))
  }))
}))
save_csv(BAL, "03_composition_by_batch_status.csv")
print(BAL, row.names=FALSE)

# Explicit comparison between original group-free residualization and
# a group-aware mean model. Use only prespecified candidate genes plus controls.
TARGET <- c("UNC13D","RAB3A","STX7","GABARAP")

# Heteroscedasticity-consistent HC3 covariance calculation (exploratory).
hc3_group <- function(fit) {
  X <- model.matrix(fit)
  j <- which(colnames(X) == "groupT2D")
  if (length(j) != 1 || qr(X)$rank < ncol(X)) return(c(beta=NA,se=NA,p=NA,lo=NA,hi=NA))
  u <- resid(fit)
  h <- hatvalues(fit)
  bread <- solve(crossprod(X))
  u2 <- (u / pmax(1-h, 1e-8))^2
  meat <- crossprod(X, X * u2)
  v <- bread %*% meat %*% bread
  se <- sqrt(max(0, v[j,j])); b <- unname(coef(fit)[j]); df <- df.residual(fit)
  if (!is.finite(se) || se == 0) return(c(beta=b,se=se,p=NA,lo=NA,hi=NA))
  c(beta=b, se=se, p=2*pt(abs(b/se),df=df,lower.tail=FALSE),
    lo=b-qt(.975,df)*se, hi=b+qt(.975,df)*se)
}

analyze_target <- function(y, data, location_group=TRUE) {
  d <- data; d$y <- as.numeric(y)
  if (location_group) {
    m <- lm(y ~ group + batch + imm + beta, data=d)
  } else {
    m <- lm(y ~ batch + imm + beta, data=d)
  }
  r <- resid(m)
  g <- d$group
  sdr <- finite_sd_ratio(r[g == "T2D"], r[g == "ND"])
  mad_nd <- mad(r[g == "ND"])
  madr <- if (is.finite(mad_nd) && mad_nd > 0) mad(r[g == "T2D"])/mad_nd else NA_real_
  z <- abs(r - ave(r, g, FUN=median))
  d$z <- z
  dm <- lm(z ~ group + batch + imm + beta, data=d)
  stat <- hc3_group(dm)
  mean_fit <- lm(y ~ group + batch + imm + beta, data=d)
  mf <- summary(mean_fit)$coefficients
  mean_coef <- unname(coef(mean_fit)["groupT2D"])
  mean_ci <- tryCatch(confint(mean_fit, "groupT2D"), error=function(e) c(NA_real_,NA_real_))
  c(sd_ratio=sdr, mad_ratio=madr,
    mean_group_effect=mean_coef, mean_ci_lo=mean_ci[1],mean_ci_hi=mean_ci[2],
    disp_group_effect=unname(stat["beta"]), disp_hc3_p=unname(stat["p"]),
    disp_hc3_ci_lo=unname(stat["lo"]), disp_hc3_ci_hi=unname(stat["hi"]))
}

cat("\n=== C. Group-aware versus group-free location: targeted sensitivity ===\n")
CAND <- do.call(rbind, lapply(TARGET, function(gene) {
  do.call(rbind, lapply(c(FALSE,TRUE), function(aware) {
    x <- analyze_target(Y[,gene], pd2, aware)
    data.frame(gene=gene, group_in_location=aware, n_ND=sum(pd2$group=="ND"),
               n_T2D=sum(pd2$group=="T2D"), t(x), row.names=NULL, check.names=FALSE)
  }))
}))
save_csv(CAND, "04_candidate_mean_dispersion_sensitivity.csv")
print(CAND, row.names=FALSE, digits=4)
cat("HC3 spread-model P values are APPROXIMATE sensitivity diagnostics, not validated confirmatory evidence.\n")

# Donor/bootstrap refits preserve joint dependence among all 45 genes.
# They stratify by disease AND batch, preserving the observed design;
# percentile CIs can be unreliable with small strata and are exploratory.
cat("\n=== D. Joint donor-bootstrap for set-level effect contrast ===\n")
STRT <- split(seq_len(nrow(pd2)), interaction(pd2$group,pd2$batch,drop=TRUE))
mean_model <- function(d) model.matrix(~ group + batch + imm + beta, data=d)
calc_set <- function(index) {
  d <- droplevels(pd2[index, , drop=FALSE]); d$group <- factor(as.character(d$group),levels=c("ND","T2D"))
  X <- mean_model(d)
  if (qr(X)$rank < ncol(X)) return(c(crino=NA,comparator=NA,difference=NA))
  R <- lm.fit(X, Y[index, , drop=FALSE])$residuals
  gr <- d$group
  ratios <- vapply(seq_len(ncol(R)), function(j) finite_sd_ratio(R[gr=="T2D",j],R[gr=="ND",j]), numeric(1))
  names(ratios) <- colnames(Y)
  a <- median(ratios[S1$gene[S1$set=="crinophagy"]],na.rm=FALSE)
  b <- median(ratios[S1$gene[S1$set=="comparator"]],na.rm=FALSE)
  c(crino=a, comparator=b, difference=a-b)
}
OBS <- calc_set(seq_len(nrow(pd2)))
if (any(!is.finite(OBS))) stop("Observed set estimates are invalid")
BOOT <- t(replicate(NBOOT, {
  ii <- unlist(lapply(STRT, function(ids) sample(ids,length(ids),replace=TRUE)),use.names=FALSE)
  calc_set(ii)
}))
valid <- apply(BOOT,1,function(z) all(is.finite(z)))
if (mean(valid) < 0.90) warning("Many bootstrap refits failed; intervals may be unstable")
BOOT <- BOOT[valid,,drop=FALSE]
if (nrow(BOOT)<100L) stop("Too few valid bootstrap replicates")
bs <- do.call(rbind,lapply(colnames(BOOT), function(v) {
  ci <- quantile(BOOT[,v],c(.025,.975),na.rm=TRUE,names=FALSE)
  data.frame(statistic=v,observed=unname(OBS[v]),ci_lo=ci[1],ci_hi=ci[2],
             n_boot=nrow(BOOT),n_failed=NBOOT-nrow(BOOT))
}))
save_csv(bs, "05_set_joint_bootstrap_CI.csv")
save_csv(as.data.frame(BOOT), "06_set_joint_bootstrap_draws.csv")
print(bs,row.names=FALSE,digits=4)
cat("CIs are conditional, exploratory donor-bootstrap intervals. They are NOT label-permutation P values.\n")

cat("\n=== E. Leave-one-donor-out influence for the four named genes ===\n")
LOO <- do.call(rbind,lapply(TARGET,function(gene) {
  full <- analyze_target(Y[,gene], pd2, TRUE)
  do.call(rbind,lapply(seq_len(nrow(pd2)),function(i) {
    ids <- setdiff(seq_len(nrow(pd2)),i)
    d <- droplevels(pd2[ids,,drop=FALSE]); d$group <- factor(as.character(d$group),levels=c("ND","T2D"))
    val <- analyze_target(Y[ids,gene],d,TRUE)
    data.frame(gene=gene,removed_donor=SAMPLES[i],removed_group=as.character(pd2$group[i]),
               sd_ratio=unname(val["sd_ratio"]),
               ratio_shift=unname(val["sd_ratio"]-full["sd_ratio"]),
               disp_hc3_p=unname(val["disp_hc3_p"]))
  }))
}))
save_csv(LOO,"07_candidate_leave_one_out.csv")
cat("Wrote ",nrow(LOO)," donor-by-gene sensitivity rows.\n")

capture.output(sessionInfo(), file=file.path(OUTPUT,"08_R_sessionInfo.txt"))
cat("\nCompleted. Files are in: ",OUTPUT,"\n",sep="")
cat("NOT addressed by this script: discovery export processing, pathology, true paired-omics,\n")
cat("candidate-selection timestamps, insulin-content links, functional mechanism and manuscript wording.\n")
