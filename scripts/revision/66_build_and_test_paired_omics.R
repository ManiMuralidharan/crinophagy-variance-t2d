# 66_build_and_test_paired_omics.R
# Referee 2: construct EXACT same-donor RNA/proteomics comparison from
# processed discovery matrices. Never insert independent-replication donors.
# Writes donor-level paired values and covariate-adjusted SD ratios with
# exploratory donor-bootstrap confidence intervals.
# RUN: source("66_build_and_test_paired_omics.R")

ROOT <- "."
OUT <- file.path(ROOT, "reviewer_validation_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
set.seed(20261113)
B <- as.integer(Sys.getenv("PAIRED_NBOOT", "1000"))
if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl")
if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
    !requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Install org.Hs.eg.db and AnnotationDbi")
DB <- org.Hs.eg.db::org.Hs.eg.db

# Recognize wide matrix: first column gene identifier, donor IDs in columns.
load_wide <- function(path, ids = c("gene_id", "feature_id", "ensembl", "ensembl_id",
                                    "gene", "gene_name", "ENTREZID", "SYMBOL")) {
  if (!file.exists(path)) stop("Missing input: ", path)
  obj <- readRDS(path)
  cat("File: ", path, "; class: ", paste(class(obj), collapse = "/"), "\n", sep = "")
  if (is.data.frame(obj)) {
    idcol <- intersect(ids, names(obj))
    if (length(idcol) > 1L) {
      # A file with multiple annotations is not safely handled by taking all
      # non-ID columns as donor measurements. Use a donor-ID signature below.
      idcol <- idcol[1L]
    }
    if (length(idcol)) {
      rid <- as.character(obj[[idcol]])
      v <- obj[, setdiff(names(obj), ids), drop = FALSE]
      x <- as.matrix(v)
      rownames(x) <- rid
    } else {
      x <- as.matrix(obj)
    }
  } else if (is.matrix(obj)) {
    x <- obj
  } else stop("Unsupported RDS class: ", paste(class(obj), collapse = ", "))
  if (is.null(rownames(x)) || is.null(colnames(x))) stop("Missing gene/donor identifiers: ", path)
  if (anyDuplicated(rownames(x))) warning("Duplicate feature IDs in ", path,
                                      "; target-level checks will require unique matches")
  suppressWarnings(storage.mode(x) <- "numeric")
  x
}

P <- load_wide(file.path(ROOT, "proteomics", "proc_prot.rds"))
R <- load_wide(file.path(ROOT, "Bulk RNA", "proc_rnaseq.rds"))
D <- read.csv(file.path(ROOT, "donor.csv"), check.names = FALSE, stringsAsFactors = FALSE)
S1 <- as.data.frame(readxl::read_excel(file.path(ROOT, "ijms-4686697-supplementary",
                                        "Supplementary_Tables.xlsx"),
                                     sheet = "S1_pre-registered_sets", skip = 2))
stopifnot(all(c("record_id", "diagnosis_computed") %in% names(D)))
if (anyDuplicated(D$record_id)) stop("Duplicate metadata donor IDs")
cat("Protein matrix:", nrow(P), "x", ncol(P), "; RNA matrix:", nrow(R), "x", ncol(R), "\n")
cat("First RNA identifiers:", paste(head(rownames(R)), collapse = ", "), "\n")
cat("First RNA sample IDs:", paste(head(colnames(R)), collapse = ", "), "\n")

# Detect transposed RNA matrices ONLY when donor metadata provides proof.
if (sum(colnames(R) %in% D$record_id) < 2 &&
    sum(rownames(R) %in% D$record_id) >= 2) {
  R <- t(R)
  cat("Transposed RNA matrix based on donor-ID matches.\n")
}
if (sum(colnames(P) %in% D$record_id) < 2) stop("Protein donor IDs cannot be linked to donor.csv")
shared <- intersect(intersect(colnames(P), colnames(R)), as.character(D$record_id))
cat("Donors in BOTH processed matrices and metadata:", length(shared), "\n")
if (length(shared) < 8L) {
  cat("Protein sample preview:", paste(head(colnames(P)), collapse = ", "), "\n")
  cat("RNA sample preview:", paste(head(colnames(R)), collapse = ", "), "\n")
  stop("Insufficient matching IDs; investigate dataset-specific donor/preparation mapping before merging")
}
D <- D[match(shared, D$record_id), , drop = FALSE]
D$group <- unname(c("None" = "ND", "Pre.T2D" = "preT2D", "Type2" = "T2D")[
  as.character(D$diagnosis_computed)])
if (anyNA(D$group)) stop("Unrecognized diagnosis codes")
D <- D[D$group %in% c("ND", "T2D"), , drop = FALSE]
cat("Two-group donors before covariate filtering:\n"); print(table(D$group))

COV <- c("donorage", "donorsex", "bodymassindex",
         "predistributionculturetime", "percentieqrecoverypostculture")
if (!all(COV %in% names(D))) stop("Discovery covariate fields missing: ",
                                 paste(setdiff(COV, names(D)), collapse = ", "))
D <- D[complete.cases(D[, COV, drop = FALSE]), , drop = FALSE]
cat("Same-donor complete-case analysis groups:\n"); print(table(D$group))
if (sum(D$group == "ND") < 4 || sum(D$group == "T2D") < 4) {
  stop("Insufficient paired complete-case donors; report this transparently")
}

# Map targets with multiple gene-ID representations (Entrez, Ensembl, symbol).
ENS <- S1[S1$gene %in% c("UNC13D", "RAB3A", "STX7", "GABARAP"),
          c("gene", "ensembl"), drop = FALSE]
annotation <- AnnotationDbi::select(DB, keys = unique(as.character(ENS$ensembl)),
                                    keytype = "ENSEMBL", columns = "ENTREZID")
annotation <- annotation[!is.na(annotation$ENTREZID), , drop = FALSE]
choose_id <- function(matrix_obj, gene, ensembl, entrez) {
  idx <- unique(which(rownames(matrix_obj) %in% c(gene, ensembl, entrez) |
                      sub("\\.[0-9]+$", "", rownames(matrix_obj)) == ensembl))
  if (!length(idx)) return(NA_integer_)
  if (length(idx) > 1L) stop("More than one possible feature for ", gene,
                           "; inspect ID mapping rather than silently averaging")
  idx
}
rows <- list(); map_audit <- list()
for (i in seq_len(nrow(ENS))) {
  gene <- as.character(ENS$gene[i]); ensembl <- as.character(ENS$ensembl[i])
  ents <- unique(annotation$ENTREZID[annotation$ENSEMBL == ensembl])
  if (length(ents) != 1L) {
    map_audit[[gene]] <- data.frame(gene = gene, ensembl = ensembl, entrez = NA_character_,
                                    protein_feature = NA_character_, rna_feature = NA_character_,
                                    status = "ambiguous_or_unmapped_entrez")
    next
  }
  pi <- choose_id(P, gene, ensembl, ents)
  ri <- choose_id(R, gene, ensembl, ents)
  status <- if (is.na(pi)) "protein_missing" else if (is.na(ri)) "rna_missing" else "paired_feature_found"
  map_audit[[gene]] <- data.frame(gene = gene, ensembl = ensembl, entrez = ents,
                                  protein_feature = if (is.na(pi)) NA_character_ else rownames(P)[pi],
                                  rna_feature = if (is.na(ri)) NA_character_ else rownames(R)[ri],
                                  status = status)
  if (is.na(pi) || is.na(ri)) next
  rows[[gene]] <- data.frame(donor = D$record_id, group = D$group, gene = gene,
                            protein_value = as.numeric(P[pi, D$record_id]),
                            rna_value = as.numeric(R[ri, D$record_id]),
                            D[, COV, drop = FALSE], check.names = FALSE)
}
write.csv(do.call(rbind, map_audit), file.path(OUT, "17_paired_feature_ID_audit.csv"), row.names = FALSE)
if (!length(rows)) stop("No target genes matched both processed matrices; inspect 17_paired_feature_ID_audit.csv")
paired <- do.call(rbind, rows)
write.csv(paired, file.path(OUT, "18_paired_omics_same_donor_values_with_covariates.csv"), row.names = FALSE)
# Also write directly to the expected location for the existing Script 61;
# omit 'batch' so it cannot apply a biologically unverified shared batch.
write.csv(paired[, c("donor", "group", "gene", "protein_value", "rna_value")],
          file.path(ROOT, "paired_omics_long.csv"), row.names = FALSE)

estimate <- function(d) {
  if (sum(d$group == "ND") < 4L || sum(d$group == "T2D") < 4L) {
    return(c(protein = NA_real_, rna = NA_real_, log_diff = NA_real_))
  }
  d$donorsex <- factor(d$donorsex)
  X <- model.matrix(~ donorage + donorsex + bodymassindex +
                    predistributionculturetime + percentieqrecoverypostculture, data = d)
  if (qr(X)$rank < ncol(X)) return(c(protein = NA_real_, rna = NA_real_, log_diff = NA_real_))
  y <- as.matrix(d[, c("protein_value", "rna_value")])
  resid <- lm.fit(X, y)$residuals
  a <- d$group == "T2D"; b <- d$group == "ND"
  srp <- sd(resid[a, 1]) / sd(resid[b, 1])
  srr <- sd(resid[a, 2]) / sd(resid[b, 2])
  if (!is.finite(srp) || !is.finite(srr) || min(srp, srr) <= 0)
    return(c(protein = NA_real_, rna = NA_real_, log_diff = NA_real_))
  c(protein = srp, rna = srr, log_diff = log(srp) - log(srr))
}

summaries <- list()
for (gene in names(rows)) {
  d <- rows[[gene]]
  d <- d[is.finite(d$protein_value) & is.finite(d$rna_value), , drop = FALSE]
  obs <- estimate(d)
  if (anyNA(obs)) {
    warning("Not enough estimable same-donor data for ", gene)
    next
  }
  strata <- split(seq_len(nrow(d)), d$group)
  bt <- t(replicate(B, {
    s <- unlist(lapply(strata, function(ids) sample(ids, length(ids), replace = TRUE)), use.names = FALSE)
    estimate(d[s, , drop = FALSE])
  }))
  valid <- apply(bt, 1, function(z) all(is.finite(z)))
  bt <- bt[valid, , drop = FALSE]
  if (nrow(bt) < 100L) {
    warning("Too few valid bootstrap refits for ", gene)
    next
  }
  ci <- apply(bt, 2, quantile, probs = c(0.025, 0.975))
  summaries[[gene]] <- data.frame(gene = gene, n_ND = sum(d$group == "ND"),
    n_T2D = sum(d$group == "T2D"), n_paired = nrow(d),
    protein_SDR = unname(obs["protein"]), protein_lo = ci[1, "protein"], protein_hi = ci[2, "protein"],
    rna_SDR = unname(obs["rna"]), rna_lo = ci[1, "rna"], rna_hi = ci[2, "rna"],
    log_protein_over_rna_SDR = unname(obs["log_diff"]),
    log_diff_lo = ci[1, "log_diff"], log_diff_hi = ci[2, "log_diff"],
    n_boot_valid = nrow(bt), n_boot_failed = B - nrow(bt))
}
if (!length(summaries)) {
  stop("No valid paired adjusted estimates; donor-level file was saved for inspection")
}
res <- do.call(rbind, summaries)
write.csv(res, file.path(OUT, "19_same_donor_covariate_adjusted_SDR.csv"), row.names = FALSE)
print(res, row.names = FALSE, digits = 4)
cat("\nDONE: 17_paired_feature_ID_audit.csv, 18_paired_omics_same_donor_values_with_covariates.csv,\n")
cat("19_same_donor_covariate_adjusted_SDR.csv, and paired_omics_long.csv.\n")
cat("NOTES: Bootstrap intervals are exploratory with few T2D donors.\n")
cat("Check Bulk RNA/proc_rnaseq_README.txt and the protein pipeline before\n")
cat("asserting exact equivalence to your originally published preprocessing.\n")
cat("Script 61 can be run on the new input but uses DIFFERENT covariate adjustment.\n")
