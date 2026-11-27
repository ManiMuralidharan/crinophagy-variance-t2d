# 65_raw_proteomics_uniprot_audit.R
# Referee 2: explicitly map unprocessed UniProt features to gene IDs;
# report raw nonmissing/positive donor counts without confusing them with
# processed/imputed measurement coverage.
# RUN: source("65_raw_proteomics_uniprot_audit.R")

ROOT <- "."
OUT <- file.path(ROOT, "reviewer_validation_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("readxl", quietly = TRUE)) stop("Install readxl: install.packages('readxl')")
if (!requireNamespace("AnnotationDbi", quietly = TRUE) ||
    !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  stop("Install org.Hs.eg.db and AnnotationDbi using BiocManager::install().")
}
DB <- org.Hs.eg.db::org.Hs.eg.db
stopifnot("UNIPROT" %in% AnnotationDbi::keytypes(DB))

read_feature_matrix <- function(path, id_field) {
  if (!file.exists(path)) stop("Missing file: ", path)
  obj <- readRDS(path)
  if (is.data.frame(obj) && id_field %in% names(obj)) {
    ids <- as.character(obj[[id_field]])
    x <- as.matrix(obj[, setdiff(names(obj), id_field), drop = FALSE])
    if (anyNA(ids) || any(!nzchar(ids))) stop("Missing feature identifiers in ", path)
    rownames(x) <- ids
  } else if (is.matrix(obj) && !is.null(rownames(obj))) {
    x <- obj
  } else {
    stop("Unexpected format in ", path, "; inspect str(readRDS(path)).")
  }
  storage.mode(x) <- "numeric"
  if (is.null(colnames(x)) || anyDuplicated(colnames(x))) stop("Invalid donor column IDs in ", path)
  x
}

RAW <- read_feature_matrix(file.path(ROOT, "proteomics", "unproc_prot.rds"), "feature_id")
PRO <- read_feature_matrix(file.path(ROOT, "proteomics", "proc_prot.rds"), "gene_id")
D <- read.csv(file.path(ROOT, "donor.csv"), check.names = FALSE, stringsAsFactors = FALSE)
S1 <- as.data.frame(readxl::read_excel(
  file.path(ROOT, "ijms-4686697-supplementary", "Supplementary_Tables.xlsx"),
  sheet = "S1_pre-registered_sets", skip = 2))
S1 <- S1[S1$set %in% c("crinophagy", "comparator"), , drop = FALSE]
stopifnot(nrow(S1) == 45L, all(c("gene", "ensembl", "set") %in% names(S1)))

cat("Raw (UniProt) features x donors:", nrow(RAW), "x", ncol(RAW), "\n")
cat("Processed (Entrez) features x donors:", nrow(PRO), "x", ncol(PRO), "\n")
cat("Shared donor identifiers:", length(intersect(colnames(RAW), colnames(PRO))), "\n")

needed <- c("record_id", "diagnosis_computed", "donorage", "donorsex",
            "bodymassindex", "predistributionculturetime", "percentieqrecoverypostculture")
if (!all(needed %in% names(D))) stop("Missing donor metadata columns: ",
                                     paste(setdiff(needed, names(D)), collapse = ", "))
if (anyDuplicated(D$record_id)) stop("Duplicate donor record_id values")
# Donor metadata may contain diagnoses beyond the ND / pre-T2D / T2D
# groups used in this analysis (or missing labels). Report them explicitly;
# do NOT guess their disease group or abort on unrelated diagnoses.
cat("\nAll diagnosis_computed labels in donor.csv:\n")
print(table(D$diagnosis_computed, useNA = "ifany"))
labels <- trimws(as.character(D$diagnosis_computed))
diagnosis_map <- c("None" = "ND", "Pre.T2D" = "preT2D", "Type2" = "T2D")
D$group <- unname(diagnosis_map[labels])
COV <- c("donorage", "donorsex", "bodymassindex", "predistributionculturetime",
         "percentieqrecoverypostculture")
shared <- D$record_id %in% intersect(colnames(RAW), colnames(PRO))
complete_cov <- complete.cases(D[, COV, drop = FALSE])
# Audit all unmatched labels among shared donors even if they are excluded.
unmatched <- shared & is.na(D$group)
if (any(unmatched)) {
  cat("\nUnmapped diagnosis labels among shared donors (NOT recoded):\n")
  print(table(D$diagnosis_computed[unmatched], useNA = "ifany"))
  write.csv(
    data.frame(record_id = D$record_id[unmatched],
               diagnosis_computed = D$diagnosis_computed[unmatched],
               reason = "unmapped_diagnosis_excluded"),
    file.path(OUT, "14_excluded_diagnosis_labels.csv"), row.names = FALSE)
  cat("Review 14_excluded_diagnosis_labels.csv to verify these exclusions.\n")
}
# Explicitly exclude unmapped diagnoses and pre-T2D donors from the
# 2-group comparison; preserve the actual ND and T2D definitions.
keep <- shared & complete_cov & !is.na(D$group) &
        D$group %in% c("ND", "T2D")
D <- D[keep, , drop = FALSE]
cat("Covariate-complete ND/T2D donors in both files:\n")
group_counts <- table(factor(D$group, levels = c("ND", "T2D")))
print(group_counts)
if (nrow(D) < 30 || !all(c("ND", "T2D") %in% D$group))
  stop("Too few ND/T2D donors; inspect metadata and excluded diagnoses")
# The submitted discovery analysis reports 87 ND and 16 T2D after
# covariate filtering. Flag any discrepancy rather than silently claiming
# its coverage counts are comparable with the manuscript.
if (unname(group_counts["ND"]) != 87L || unname(group_counts["T2D"]) != 16L)
  stop("Discovery group counts differ from the manuscript's expected 87 ND / 16 T2D; ",
       "review the diagnosis exclusion file and covariate filtering before proceeding")

# Source gene annotation: Ensembl -> Entrez, never silently choose among multiple IDs.
valid_ens <- intersect(unique(as.character(S1$ensembl)),
                       AnnotationDbi::keys(DB, keytype = "ENSEMBL"))
ens_map <- unique(AnnotationDbi::select(DB, keys = valid_ens,
                                        keytype = "ENSEMBL", columns = "ENTREZID"))
ens_map <- ens_map[!is.na(ens_map$ENTREZID), , drop = FALSE]

# Raw protein identifiers are UniProt accessions (e.g. A0A024RBG1).
# A suffix such as '-2' is a protein isoform indicator: preserve original,
# remove suffix only for annotation lookup. Unmapped features remain unmapped.
raw_original <- rownames(RAW)
raw_lookup <- sub("-[0-9]+$", "", raw_original)
known_uniprot <- AnnotationDbi::keys(DB, keytype = "UNIPROT")
valid_raw <- intersect(unique(raw_lookup), known_uniprot)
if (!length(valid_raw)) stop("No UniProt accessions recognized; verify identifier format")
raw_map <- unique(AnnotationDbi::select(DB, keys = valid_raw,
                                        keytype = "UNIPROT", columns = "ENTREZID"))
raw_map <- raw_map[!is.na(raw_map$ENTREZID), , drop = FALSE]
raw_keys_to_entrez <- split(raw_map$ENTREZID, raw_map$UNIPROT)
raw_keys_to_entrez <- lapply(raw_keys_to_entrez, unique)
raw_unique <- vapply(raw_lookup, function(k) {
  ids <- raw_keys_to_entrez[[k]]
  if (is.null(ids) || length(ids) != 1L) NA_character_ else ids
}, character(1))
raw_mapping_status <- vapply(raw_lookup, function(k) {
  ids <- raw_keys_to_entrez[[k]]
  if (is.null(ids)) "UNMAPPED" else if (length(ids) == 1L) "UNIQUE" else "AMBIGUOUS"
}, character(1))
write.csv(data.frame(feature_id = raw_original, lookup_accession = raw_lookup,
                     mapped_entrez = raw_unique, status = raw_mapping_status),
          file.path(OUT, "15_raw_uniprot_entrez_mapping.csv"), row.names = FALSE)
cat("Raw UniProt mappings: ", paste(names(table(raw_mapping_status)),
                                    as.integer(table(raw_mapping_status)), collapse = "; "), "\n", sep = "")

panel_results <- lapply(seq_len(nrow(S1)), function(i) {
  ent <- unique(ens_map$ENTREZID[ens_map$ENSEMBL == as.character(S1$ensembl[i])])
  raw_rows <- if (length(ent) == 1L) which(raw_unique == ent) else integer(0)
  pro_rows <- if (length(ent) == 1L) which(rownames(PRO) == ent) else integer(0)
  n_nd <- sum(D$group == "ND"); n_t2d <- sum(D$group == "T2D")
  count_nonmiss <- function(X, ix, ids, positive = FALSE) {
    if (!length(ix)) return(NA_integer_)
    v <- X[ix, ids, drop = FALSE]
    flag <- if (positive) is.finite(v) & v > 0 else is.finite(v)
    sum(colSums(flag) > 0)
  }
  nd <- D$record_id[D$group == "ND"]
  t2 <- D$record_id[D$group == "T2D"]
  data.frame(gene = as.character(S1$gene[i]), set = as.character(S1$set[i]),
             ensembl = as.character(S1$ensembl[i]),
             entrez = if (length(ent) == 1L) ent else NA_character_,
             annotation_status = if (length(ent) == 1L) "unique_ensembl_entrez"
                                 else if (!length(ent)) "unmapped_ensembl" else "ambiguous_ensembl",
             raw_accessions = length(raw_rows), processed_rows = length(pro_rows),
             eligible_ND = n_nd, eligible_T2D = n_t2d,
             raw_nonmissing_ND = count_nonmiss(RAW, raw_rows, nd),
             raw_nonmissing_T2D = count_nonmiss(RAW, raw_rows, t2),
             raw_positive_ND = count_nonmiss(RAW, raw_rows, nd, TRUE),
             raw_positive_T2D = count_nonmiss(RAW, raw_rows, t2, TRUE),
             processed_nonmissing_ND = count_nonmiss(PRO, pro_rows, nd),
             processed_nonmissing_T2D = count_nonmiss(PRO, pro_rows, t2))
})
res <- do.call(rbind, panel_results)
write.csv(res, file.path(OUT, "16_raw_and_processed_panel_coverage.csv"), row.names = FALSE)
print(res[res$gene %in% c("UNC13D", "RAB3A", "STX7", "GABARAP"), ],
      row.names = FALSE)
cat("\nDONE: 15_raw_uniprot_entrez_mapping.csv and 16_raw_and_processed_panel_coverage.csv\n")
cat("CAUTION: raw nonmissing is a feature-availability count, NOT proof of\n",
    "biological presence. Multi-accession proteins were combined only for\n",
    "donor nonmissing counts. Ambiguous UniProt mappings were excluded.\n",
    "Processed finite counts can include imputation.\n", sep = "")
