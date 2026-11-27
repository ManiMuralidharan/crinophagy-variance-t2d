# 62_discovery_coverage_audit.R
# Referee 2: discovery source/protein IDs, detection, n by disease group,
# covariate completeness, exploratory variance estimates in processed data.
# Does not by itself verify a preprint's preprocessing pipeline, prove
# imputation happened, or establish the original donor reclassification.

ROOT <- "."
OUT <- file.path(ROOT,"reviewer_validation_output")
XLS <- file.path(ROOT,"ijms-4686697-supplementary","Supplementary_Tables.xlsx")
PFILE <- file.path(ROOT,"proteomics","proc_prot.rds")
UFILE <- file.path(ROOT,"proteomics","unproc_prot.rds")
DFILE <- file.path(ROOT,"donor.csv")
for (p in c(XLS,PFILE,DFILE))if(!file.exists(p))stop("Missing required input: ",p)
if(!requireNamespace("readxl",quietly=TRUE))stop("install.packages('readxl')")
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)

s1_all <- as.data.frame(readxl::read_excel(XLS,sheet="S1_pre-registered_sets",skip=2))
if (!all(c("gene","ensembl","set") %in% names(s1_all)))stop("Unexpected S1 columns")
# Exclude STX1A (post-hoc) from the 45 frozen gene-set members.
s1 <- s1_all[!is.na(s1_all$set) &
             s1_all$set %in% c("crinophagy", "comparator"), , drop=FALSE]
stopifnot(nrow(s1) == 45L, sum(s1$set == "crinophagy") == 31L,
          sum(s1$set == "comparator") == 14L)
cat("Panel membership: 31 crinophagy + 14 comparator; excluded",
    nrow(s1_all)-nrow(s1), "non-frozen row(s)\n")
don <- read.csv(DFILE,stringsAsFactors=FALSE,check.names=FALSE)
needed <- c("record_id","diagnosis_computed","donorage","donorsex",
            "bodymassindex","predistributionculturetime","percentieqrecoverypostculture")
if (!all(needed %in% names(don)))stop("Discovery donor metadata has missing expected columns: ",
                                      paste(setdiff(needed,names(don)),collapse=", "))

to_matrix <- function(path) {
  obj <- readRDS(path)
  if (is.data.frame(obj) && "gene_id" %in% names(obj)) {
    ids <- as.character(obj$gene_id)
    x <- as.matrix(obj[,setdiff(names(obj),"gene_id"),drop=FALSE])
    rownames(x) <- ids
  } else if (is.matrix(obj) || is.data.frame(obj)) {
    x <- as.matrix(obj)
  } else stop("Unsupported .rds object structure in ",path)
  if(is.null(rownames(x)) || is.null(colnames(x)))stop("Matrix lacks protein IDs or donor IDs: ",path)
  if(anyDuplicated(rownames(x)))stop("Repeated protein identifiers in ",path,"; resolve aggregation explicitly")
  storage.mode(x) <- "numeric"
  x
}
P <- to_matrix(PFILE)
# Unprocessed exports may have a different structure (long tables, nested lists,
# or identifiers stored in columns). Do NOT infer their orientation or IDs.
# Keep the processed-data audit running; record the raw-data blocker separately.
raw_status <- "missing_unprocessed_file"
U <- NULL
if (file.exists(UFILE)) {
  U <- tryCatch(to_matrix(UFILE), error=function(e) {
    raw_status <<- paste("unprocessed_structure_unresolved:", conditionMessage(e))
    warning("Raw-protein coverage NOT evaluated: ", conditionMessage(e))
    # Structural diagnostics only; no expression values or inferred identifiers.
    raw_obj <- tryCatch(readRDS(UFILE), error=function(err) NULL)
    structure_text <- if (is.null(raw_obj)) {
      "Unable to inspect unprocessed file"
    } else {
      c(paste("Class:", paste(class(raw_obj), collapse="/")),
        paste("Dimensions:", paste(dim(raw_obj), collapse=" x ")),
        paste("Column names:", paste(head(colnames(raw_obj), 35), collapse=" | ")),
        paste("Element names:", paste(head(names(raw_obj), 35), collapse=" | ")),
        paste("Rowname examples:", paste(head(rownames(raw_obj), 5), collapse=" | ")))
    }
    writeLines(structure_text, file.path(OUT, "10_unprocessed_structure_diagnostic.txt"))
    NULL
  })
  if (!is.null(U)) raw_status <- "unprocessed_matrix_parsed"
}
writeLines(raw_status, file.path(OUT, "10_unprocessed_coverage_status.txt"))

if (!all(colnames(P) %in% don$record_id))stop("Processed protein columns cannot all match donor.csv")
don <- don[match(colnames(P),don$record_id),,drop=FALSE]
stopifnot(identical(as.character(don$record_id),colnames(P)))
don$grp <- unname(c("None"="ND","Pre.T2D"="preT2D","Type2"="T2D")[as.character(don$diagnosis_computed)])
if(anyNA(don$grp))stop("Unmapped diagnosis_computed values; inspect original donor metadata")
COV <- c("donorage","donorsex","bodymassindex",
         "predistributionculturetime","percentieqrecoverypostculture")
keep <- don$grp %in% c("ND","T2D") & complete.cases(don[,COV,drop=FALSE])
d <- don[keep,,drop=FALSE]
cat("\nProcessed file:",PFILE," dimensions:",nrow(P),"x",ncol(P),"\n")
print(table(d$grp))

# Explicit accession mapping: use matching IDs only. A mapping with >1
# candidate protein in the source is labelled ambiguous, never silently picked.
map <- data.frame(gene=as.character(s1$gene),ensembl=as.character(s1$ensembl),
                  entrez=NA_character_, mapping_status="not_mapped")
if(requireNamespace("org.Hs.eg.db",quietly=TRUE) &&
   requireNamespace("AnnotationDbi",quietly=TRUE)) {
  ref <- AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                               keys=unique(map$ensembl),keytype="ENSEMBL",
                               columns="ENTREZID")
  for(i in seq_len(nrow(map))) {
    cand <- unique(as.character(na.omit(ref$ENTREZID[ref$ENSEMBL==map$ensembl[i]])))
    hits <- intersect(cand,rownames(P))
    if(length(hits)==1L) {
      map$entrez[i] <- hits
      map$mapping_status[i] <- "one_protein_id_matched"
    } else if(length(hits)>1L) {
      map$mapping_status[i] <- "AMBIGUOUS_multiple_processed_IDs"
    } else {
      map$mapping_status[i] <- "NOT_FOUND_in_processed_IDs"
    }
  }
} else warning("Bioconductor AnnotationDbi/org.Hs.eg.db unavailable; mapping cannot be completed.")
write.csv(map,file.path(OUT,"11_discovery_gene_identifier_mapping.csv"),row.names=FALSE)

# NA means no conclusive match, NOT proof a protein was not measured.
res <- do.call(rbind,lapply(seq_len(nrow(map)), function(i) {
  hit <- map$entrez[i]
  p_nd <- p_t2d <- u_nd <- u_t2d <- NA_integer_
  p_sdr <- NA_real_
  unproc_status <- "not_evaluable_ID_or_file"
  if(!is.na(hit)) {
    y <- as.numeric(P[hit,don$record_id])
    p_nd <- sum(is.finite(y[don$grp=="ND" & keep]))
    p_t2d <- sum(is.finite(y[don$grp=="T2D" & keep]))
    yi <- as.numeric(P[hit,d$record_id])
    complete <- is.finite(yi)
    if(sum(complete & d$grp=="ND")>=4 && sum(complete & d$grp=="T2D")>=4) {
      dm <- d[complete,,drop=FALSE]
      dm$donorsex <- factor(dm$donorsex)
      X <- model.matrix(~ donorage + donorsex + bodymassindex +
                         predistributionculturetime + percentieqrecoverypostculture,data=dm)
      if(qr(X)$rank==ncol(X)) {
        rr <- lm.fit(X,yi[complete])$residuals
        n <- dm$grp=="ND"; t <- dm$grp=="T2D"
        if(sd(rr[n])>0)p_sdr <- sd(rr[t])/sd(rr[n])
      }
    }
    # Raw ID matches cannot be assumed across Entrez and UniProt exports.
    if(!is.null(U) && hit %in% rownames(U) && all(d$record_id %in% colnames(U))) {
      unproc_status <- "entrez_id_matched_raw"
      raw <- as.numeric(U[hit,d$record_id])
      u_nd <- sum(is.finite(raw[d$grp=="ND"]))
      u_t2d <- sum(is.finite(raw[d$grp=="T2D"]))
    }
  }
  data.frame(gene=map$gene[i],set=s1$set[i],ensembl=map$ensembl[i],
             processed_id=hit,mapping_status=map$mapping_status[i],
             proc_finite_ND=p_nd,proc_finite_T2D=p_t2d,
             raw_finite_ND=u_nd,raw_finite_T2D=u_t2d,
             raw_id_status=unproc_status,
             processed_SDR_same_covariates=p_sdr)
}))
write.csv(res,file.path(OUT,"12_discovery_panel_coverage.csv"),row.names=FALSE)
print(res[res$gene %in% c("UNC13D","RAB3A","GABARAP","STX7"),],row.names=FALSE)

source_paths <- c(PFILE,DFILE,XLS)
if(file.exists(UFILE))source_paths<-c(source_paths,UFILE)
files <- data.frame(path=source_paths,bytes=as.numeric(file.info(source_paths)$size),
                    md5=unname(tools::md5sum(source_paths)))
write.csv(files,file.path(OUT,"13_source_file_fingerprints.csv"),row.names=FALSE)

# Identify potential pathology fields WITHOUT inferring unavailable diagnoses.
path_cols <- grep("tumou|tumor|neoplas|malig|cancer|patholog|pancrea|surg|diagnos",
                  names(don),ignore.case=TRUE,value=TRUE)
write.csv(data.frame(column=path_cols),file.path(OUT,"14_available_pathology_columns.csv"),row.names=FALSE)
cat("\nSaved coverage, mapping, file fingerprints, and potential pathology fields.\n")
cat("Processed finite values MAY be imputed: they are not raw detection frequencies.\n")
cat("If raw identifiers or layout are unresolved, raw counts stay NA\n")
cat("until an audited UniProt-to-Entrez mapping is supplied.\n")
cat("Neither processing steps nor 21 status reclassifications can be verified from these files alone.\n")
