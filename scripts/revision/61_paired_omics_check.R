# 61_paired_omics_check.R
# Referee 2: reconcile same-donor paired discovery protein/RNA SD ratios
# with independent replication. DOES NOT reconstruct missing input data.
# REQUIRED CSV (already normalized on suitable quantitative scales):
#   paired_omics_long.csv: donor,group,gene,protein_value,rna_value
# Optional column: batch (use identical donor batch covariate in both layers).
# Each row is one donor-gene pair. Include UNC13D and RAB3A (and optionally
# STX7, GABARAP). Do not mix discovery and independent-replication donors.

ROOT <- "."
INFILE <- file.path(ROOT,"paired_omics_long.csv")
OUT <- file.path(ROOT,"reviewer_validation_output")
NBOOT <- as.integer(Sys.getenv("REVIEWER_NBOOT", "1200"))
set.seed(20261113)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

if (!file.exists(INFILE)) {
  template <- data.frame(donor=character(),group=character(),gene=character(),
                         protein_value=numeric(),rna_value=numeric(),batch=character())
  write.csv(template,file.path(OUT,"PAIRED_OMICS_INPUT_TEMPLATE.csv"),row.names=FALSE)
  stop("Cannot calculate paired-donor results without matched values. ",
       "Create paired_omics_long.csv using the template in ",OUT,
       ". Each row must have both modalities measured in the same donor.")
}
dat <- read.csv(INFILE,stringsAsFactors=FALSE,check.names=FALSE)
req <- c("donor","group","gene","protein_value","rna_value")
missing <- setdiff(req,names(dat))
if(length(missing))stop("Missing columns: ",paste(missing,collapse=", "))
dat$protein_value <- suppressWarnings(as.numeric(dat$protein_value))
dat$rna_value <- suppressWarnings(as.numeric(dat$rna_value))
dat$group <- factor(dat$group,levels=c("ND","T2D"))
if (anyNA(dat[,c("donor","gene","group")])) stop("Missing IDs/gene/status")
if (anyDuplicated(paste(dat$donor,dat$gene,sep="||"))) stop("Duplicate donor/gene rows: resolve technical replicates explicitly")
if ("batch" %in% names(dat) && anyNA(dat$batch)) stop("Missing batch values: repair or remove the column deliberately")

calc_ratio <- function(y,grp) {
  a <- y[grp=="T2D"]; b <- y[grp=="ND"]
  if (length(a)<4 || length(b)<4 || sd(b)==0) return(NA_real_)
  sd(a)/sd(b)
}
calc_both <- function(d) {
  # A group-aware location model separates differences in mean from spread.
  # Batch adjustment is possible only when its design is estimable.
  X <- if ("batch" %in% names(d) && length(unique(d$batch)) > 1)
         model.matrix(~ group + factor(batch),data=d)
       else model.matrix(~ group,data=d)
  if (qr(X)$rank < ncol(X)) return(c(protein=NA,rna=NA,log_ratio_difference=NA))
  ep <- lm.fit(X,d$protein_value)$residuals
  er <- lm.fit(X,d$rna_value)$residuals
  p <- calc_ratio(ep,d$group); r <- calc_ratio(er,d$group)
  c(protein=p,rna=r,log_ratio_difference=log(p)-log(r))
}

results <- list(); draws <- list()
for (gene in intersect(c("UNC13D","RAB3A","STX7","GABARAP"),unique(dat$gene))) {
  d <- dat[dat$gene==gene,,drop=FALSE]
  n_before <- nrow(d)
  d <- d[is.finite(d$protein_value)&is.finite(d$rna_value),,drop=FALSE]
  d$group <- droplevels(d$group)
  if (length(unique(d$group)) != 2 || any(table(d$group)<4)) {
    warning(gene,": insufficient overlapping paired data; omitted")
    next
  }
  # Exact same donors in both layers; bootstrap jointly within groups.
  obs <- calc_both(d)
  s <- split(seq_len(nrow(d)),d$group)
  B <- t(replicate(NBOOT, {
    ix <- unlist(lapply(s,function(ii) sample(ii,length(ii),replace=TRUE)),use.names=FALSE)
    calc_both(d[ix,,drop=FALSE])
  }))
  ok <- apply(B,1,function(z) all(is.finite(z)))
  B <- B[ok,,drop=FALSE]
  if(nrow(B)<100)stop("Insufficient valid bootstrap draws for ",gene)
  CI <- apply(B,2,quantile,probs=c(.025,.975),na.rm=TRUE)
  results[[gene]] <- data.frame(gene=gene,n_paired_before_missingness=n_before,
     n_paired_used=nrow(d),n_ND=sum(d$group=="ND"),n_T2D=sum(d$group=="T2D"),
     protein_sd_ratio=obs["protein"],protein_ci_lo=CI[1,"protein"],protein_ci_hi=CI[2,"protein"],
     rna_sd_ratio=obs["rna"],rna_ci_lo=CI[1,"rna"],rna_ci_hi=CI[2,"rna"],
     log_protein_to_rna_ratio=obs["log_ratio_difference"],
     log_diff_ci_lo=CI[1,"log_ratio_difference"],log_diff_ci_hi=CI[2,"log_ratio_difference"],
     n_boot=nrow(B))
  draws[[gene]] <- data.frame(gene=gene,B,check.names=FALSE)
}
if(!length(results))stop("No genes could be analyzed; check paired input")
out <- do.call(rbind,results)
write.csv(out,file.path(OUT,"09_same_donors_paired_omics.csv"),row.names=FALSE)
write.csv(do.call(rbind,draws),file.path(OUT,"10_paired_omics_bootstrap_draws.csv"),row.names=FALSE)
print(out,row.names=FALSE,digits=4)
cat("\nInterpretation: compare these SAME-DONOR figures to discovery proteomics\n")
cat("(published 2.065 for UNC13D, 1.725 for RAB3A) AND INDEPENDENT replication\n")
cat("(published 2.001, 1.101), in separate columns of the final manuscript table.\n")
cat("CI methods are percentile bootstrap and are exploratory given small T2D cells.\n")
cat("This script does not validate original normalization, cell source or dataset independence.\n")
