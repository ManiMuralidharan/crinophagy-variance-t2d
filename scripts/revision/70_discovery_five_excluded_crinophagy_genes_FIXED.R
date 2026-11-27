# =====================================================================
# 70_discovery_five_excluded_crinophagy_genes.R
# READ-ONLY AUDIT of the 31 frozen crinophagy genes vs 26 discovery tests.
# Requires outputs from Script 65 and historical per_protein_variance.csv.
# IMPORTANT: panel absence and raw missingness do not establish the exact
#             provider filtering rule (release-specific source unverified).
# =====================================================================
ROOT <- "."
OUT <- file.path(ROOT, "reviewer_validation_output")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
if (!requireNamespace("readxl", quietly = TRUE)) stop("Run install.packages('readxl')")

choose_one <- function(paths, label) {
  paths <- unique(paths[file.exists(paths)])
  if (!length(paths)) stop(label, " not found. Checked:\n", paste(paths, collapse = "\n"))
  paths[1]
}
xl <- choose_one(c(file.path(ROOT, "ijms-4686697-supplementary", "Supplementary_Tables.xlsx"),
                   file.path(ROOT, "Supplementary_Tables.xlsx")), "Supplementary_Tables.xlsx")
pp_path <- file.path(ROOT, "per_protein_variance.csv")
cov_path <- file.path(OUT, "16_raw_and_processed_panel_coverage.csv")
map_path <- file.path(OUT, "15_raw_uniprot_entrez_mapping.csv")
raw_path <- file.path(ROOT, "proteomics", "unproc_prot.rds")
for (f in c(pp_path, cov_path, map_path, raw_path)) if (!file.exists(f)) stop("Missing file: ", f)

s1 <- as.data.frame(readxl::read_excel(xl, sheet = "S1_pre-registered_sets", skip = 2))
s1$gene <- trimws(as.character(s1$gene)); s1$set <- trimws(as.character(s1$set))
ref <- s1[s1$set == "crinophagy", c("gene", "ensembl"), drop = FALSE]
if (nrow(ref) != 31 || anyDuplicated(ref$gene)) stop("Expected 31 distinct crinophagy symbols from frozen Table S1; inspect file")
pp <- read.csv(pp_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("symbol", "set", "sd_ratio", "p", "q") %in% names(pp))) stop("Unexpected per_protein_variance.csv columns")
tested <- unique(trimws(as.character(pp$symbol[pp$set == "crinophagy"])))
if (length(tested) != 26) stop("Expected 26 tested crinophagy proteins; found ",length(tested))
if (!all(tested %in% ref$gene)) stop("Tested symbols not in frozen S1: ", paste(setdiff(tested,ref$gene),collapse=", "))
missing <- setdiff(ref$gene, tested)
if (length(missing) != 5) stop("Expected exactly five excluded genes; found ",length(missing))

cov <- read.csv(cov_path, stringsAsFactors = FALSE, check.names = FALSE)
needed <- c("gene","set","annotation_status","entrez","raw_accessions","processed_rows",
            "raw_nonmissing_ND","raw_nonmissing_T2D","eligible_ND","eligible_T2D")
if (!all(needed %in% names(cov))) stop("Missing coverage fields: ", paste(setdiff(needed, names(cov)),collapse=", "))
cov <- cov[cov$set == "crinophagy", ,drop=FALSE]
if (anyDuplicated(cov$gene)) stop("Duplicate crinophagy symbols in coverage output")
idx <- match(ref$gene,cov$gene)
if (anyNA(idx)) stop("Coverage output missing: ",paste(ref$gene[is.na(idx)],collapse=", "))
audit <- cbind(ref, cov[idx,setdiff(names(cov),c("gene","set","ensembl")),drop=FALSE])
audit$in_historical_discovery_test <- audit$gene %in% tested
if (any(audit$in_historical_discovery_test & audit$processed_rows == 0, na.rm=TRUE)) {
  warning("A tested gene is reported absent from processed matrix. Check mapping before use.")
}
audit$observation <- ifelse(
  audit$in_historical_discovery_test, "TESTED in original processed-matrix output",
  ifelse(is.na(audit$entrez) | audit$annotation_status != "unique_ensembl_entrez",
         "NOT TESTED; Ensembl-to-Entrez mapping unresolved",
         ifelse(is.na(audit$processed_rows), "NOT TESTED; processed mapping unresolved",
                ifelse(audit$processed_rows == 0,
                       "NOT TESTED; gene absent from processed matrix",
                       "NOT TESTED despite processed-row match; inspect original code"))))
audit$raw_detection_rate_ND <- with(audit, raw_nonmissing_ND / eligible_ND)
audit$raw_detection_rate_T2D <- with(audit, raw_nonmissing_T2D / eligible_T2D)
audit$exact_provider_filter_reason_verified <- FALSE
write.csv(audit, file.path(OUT,"38_crinophagy_discovery_31_gene_coverage.csv"),row.names=FALSE)

# Inspect raw accession-level NA fractions; use all 129 raw donors because
# source-level filtering may not have used ND/T2D donor subsets.
map <- read.csv(map_path,stringsAsFactors=FALSE,check.names=FALSE)
stopifnot(all(c("feature_id","mapped_entrez","status") %in% names(map)))
raw_obj <- readRDS(raw_path)
if (!is.data.frame(raw_obj) || !"feature_id" %in% names(raw_obj)) stop("Unexpected unproc_prot.rds layout")
raw_id <- as.character(raw_obj$feature_id)
if (anyDuplicated(raw_id)) stop("Duplicate raw feature IDs: inspect")
rx <- as.matrix(raw_obj[,setdiff(names(raw_obj),"feature_id"), drop=FALSE])
storage.mode(rx) <- "numeric"
rownames(rx) <- raw_id
five_indices <- which(!audit$in_historical_discovery_test)
accession_details <- lapply(five_indices, function(i) {
  z <- audit[i,]
  ent <- as.character(z$entrez)
  if (is.na(z$entrez) || !nzchar(ent)) return(data.frame(gene=z$gene, entrez=NA_character_,
    raw_accession=NA_character_, raw_na_fraction_all_donors=NA_real_,
    raw_nonmissing_all_donors=NA_integer_, note="No unique Entrez ID; cannot identify raw accession"))
  acc <- unique(as.character(map$feature_id[!is.na(map$mapped_entrez) & map$mapped_entrez==ent & map$status=="UNIQUE"]))
  if (!length(acc)) return(data.frame(gene=z$gene,entrez=ent,raw_accession=NA_character_,
    raw_na_fraction_all_donors=NA_real_,raw_nonmissing_all_donors=NA_integer_,
    note="No uniquely mapped raw UniProt accession; may be annotation issue"))
  lapply(acc,function(id) {
    if (!id %in% rownames(rx)) stop("Raw accession missing from matrix: ",id)
    v <- rx[id,]
    data.frame(gene=z$gene,entrez=ent,raw_accession=id,
      raw_na_fraction_all_donors=mean(!is.finite(v)),
      raw_nonmissing_all_donors=sum(is.finite(v)),
      note="Raw nonmissing counts do not establish biological presence or exclusion mechanism")
  })
})
# Flatten only the five excluded-gene results. Every returned element is
# either a single data.frame (no accession) or a list of data.frames
# (one per accession). Never pass NULL to do.call().
accession_rows <- list()
for (z in accession_details) {
  if (is.data.frame(z)) {
    accession_rows[[length(accession_rows) + 1L]] <- z
  } else if (is.list(z) && all(vapply(z, is.data.frame, logical(1)))) {
    for (item in z) accession_rows[[length(accession_rows) + 1L]] <- item
  } else {
    stop("Unexpected accession-result type; inspect raw identifier mapping")
  }
}
if (!length(accession_rows)) stop("No accession summary rows generated")
accession_details <- do.call(rbind, accession_rows)
rownames(accession_details) <- NULL
write.csv(accession_details,file.path(OUT,"39_five_excluded_raw_accessions.csv"),row.names=FALSE)

five <- audit[!audit$in_historical_discovery_test,,drop=FALSE]
cat("\nFIVE EXCLUDED GENES (against historical processed-matrix test output):\n")
print(five[,c("gene","entrez","annotation_status","raw_accessions","processed_rows",
               "raw_nonmissing_ND","raw_nonmissing_T2D","observation")],row.names=FALSE)
cat("\nRAW ACCESSION DETAILS:\n")
print(accession_details,row.names=FALSE)
cat("\nPotential 50% missingness rule: inspect raw_na_fraction_all_donors.\n")
cat("A feature having <=50% raw missingness but absent from processed export\n")
cat("means simple >50%-missing filtering alone cannot explain its absence.\n")
cat("Even if >50%, coincidence does not prove exact original filtering.\n")
summary_txt <- c(
 paste0("Reference set: ",length(ref)," genes; historical processed-matrix tests: ",length(tested),"; excluded: ",length(missing)),
 paste0("Excluded symbols: ",paste(sort(missing),collapse=", ")),
 "Per-gene statuses are observable facts, not verified source-processing causes.",
 "Absent from processed matrix does NOT equal unmeasured in raw proteomics.",
 "Raw NA fractions are calculated over all 129 raw donor columns; zeros are not recoded as NA.",
 "Exact stage of provider filtering requires matching release-specific processing code.",
 "Original discovery P/q values are exploratory until residual-test calibration and model assumptions are addressed.")
writeLines(summary_txt,file.path(OUT,"40_excluded_crinophagy_interpretation.txt"))
cat("\nWrote 38_crinophagy_discovery_31_gene_coverage.csv,\n")
cat("39_five_excluded_raw_accessions.csv, and 40_excluded_crinophagy_interpretation.txt in ",OUT,"\n",sep="")
