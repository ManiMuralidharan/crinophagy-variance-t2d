# 63_bootstrap_coverage_sensitivity.R
# Referee 2: approximate interval coverage at the actual ND/T2D group sizes.
# Compares percentile bootstrap with analytical log-SD-ratio CI, under
# NORMAL and HEAVY-TAILED simulation models with known SD ratios.
# These simulations DO NOT include nuisance regression, selection, missing
# protein measurements, imputation, batch, or beta/immune confounding.
# Thus they cannot validate the published adjusted intervals on their own.

ROOT <- "."
OUT <- file.path(ROOT,"reviewer_validation_output")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
set.seed(20261113)
NSIM <- as.integer(Sys.getenv("REVIEWER_NSIM", "350"))
NBOOT <- as.integer(Sys.getenv("REVIEWER_COVER_BOOT", "300"))
if (NSIM < 50 || NBOOT < 100) stop("Need at least 50 simulations and 100 bootstraps")

sim_sample <- function(n, scale, kind) {
  x <- if (kind=="normal") rnorm(n) else rt(n, df=5)/sqrt(5/3)
  scale*x
}
ratio <- function(x,y) sd(x)/sd(y)
run_one <- function(n1,n0,true,dist) {
  covered_pct <- covered_analytical <- logical(NSIM)
  width_pct <- numeric(NSIM)
  for(i in seq_len(NSIM)) {
    a <- sim_sample(n1,true,dist); b <- sim_sample(n0,1,dist)
    t <- ratio(a,b)
    bt <- replicate(NBOOT,ratio(sample(a,n1,TRUE),sample(b,n0,TRUE)))
    ci <- unname(quantile(bt,c(.025,.975),na.rm=TRUE))
    # Approximate log-ratio SE is accurate chiefly under normality.
    se <- sqrt(1/(2*(n1-1)) + 1/(2*(n0-1)))
    ca <- exp(log(t) + c(-1,1)*qnorm(.975)*se)
    covered_pct[i] <- ci[1]<=true && true<=ci[2]
    covered_analytical[i] <- ca[1]<=true && true<=ca[2]
    width_pct[i] <- diff(ci)
  }
  data.frame(n_T2D=n1,n_ND=n0,true_ratio=true,distribution=dist,
             percentile_coverage=mean(covered_pct),
             normal_logratio_coverage=mean(covered_analytical),
             median_percentile_width=median(width_pct),
             n_sim=NSIM,n_boot=NBOOT)
}
DESIGNS <- list(replication=c(36,18), discovery_maximum=c(16,87))
out <- do.call(rbind,lapply(names(DESIGNS),function(study) {
  do.call(rbind,lapply(c("normal","t5"),function(d) {
    do.call(rbind,lapply(c(1,1.5,2),function(ra) {
      z <- run_one(DESIGNS[[study]][1],DESIGNS[[study]][2],ra,d)
      z$study <- study; z
    }))
  }))
}))
out <- out[,c("study","n_T2D","n_ND","true_ratio","distribution",
             "percentile_coverage","normal_logratio_coverage",
             "median_percentile_width","n_sim","n_boot")]
write.csv(out,file.path(OUT,"15_basic_SD_ratio_interval_coverage.csv"),row.names=FALSE)
print(out,row.names=FALSE,digits=3)
cat("\nThese calibrate unadjusted SD-ratio intervals under specified distributions.\n")
cat("They do NOT establish coverage for the manuscript's covariate-adjusted bootstrap.\n")
