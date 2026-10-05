# Inter-individual variance of crinophagy machinery in type 2 diabetic human islets

Analysis code for the study of inter-individual variance, rather than
mean abundance, in crinophagy-related proteins in human pancreatic
islets from donors with and without type 2 diabetes.

**Mani M, Martin TFJ.** Department of Biochemistry, University of
Wisconsin-Madison.

---

## What this repository contains

Scripts only. No data files are included: the primary data belong to
HumanIslets.com and to the NCBI Gene Expression Omnibus, and must be
obtained from those sources (see **Data** below).

```
genesets/      Pre-registered gene sets, with version history
discovery/     Discovery cohort analysis (HumanIslets proteomics)
replication/   Replication cohort analysis (GSE164416 RNA-seq)
specificity/   Comparator and specificity panels
figures/       Figure generation
tables/        Supplementary table generation
```

Scripts are numbered in run order. The numbering is continuous across
folders, so `08` in `replication/` follows `07` in `discovery/`.

---

## Pre-registered gene sets

Two protein sets were defined before any statistical evaluation and
deposited with a timestamped DOI:

**https://doi.org/10.5281/zenodo.23060343**

- `frozen_gene_sets.txt` — version 1
- `frozen_gene_sets_v2_clean.txt` — version 2, used for all analyses
- `frozen_gene_sets_v2_header.txt` — record of the single revision,
  with its date and reason

The comparator set was revised once before any statistical test:
MAP1LC3B and ATG12 were not detected in the discovery proteomics, so
the set was expanded to the full ATG8 family and the selective
autophagy cargo receptors. ATG12 was excluded as it was detected in
neither dataset. The header file documents this.

---

## Data

Neither dataset is redistributed here.

**Discovery cohort.** Islet proteomics, bulk RNA-seq and donor
metadata from [HumanIslets.com](https://www.humanislets.com).
Registration is required. Download the proteomics (`proc_prot.rds`,
`unproc_prot.rds`), bulk RNA-seq (`proc_rnaseq.rds`) and the donor
list (`donor.csv`).

**Replication cohort.** GEO accession
[GSE164416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164416)
(Wigger et al., *Nat Metab* 2021). Retrieved with GEOquery; see
`replication/08_gse164416_inspect.py`.

**Single-cell reference.** Mean-difference values were taken from the
published supplementary tables of Lawlor et al., *Genome Res* 2017
(GSE86473). No reanalysis of the primary single-cell data was
performed.

---

## Environment

- **R 4.6.1** (2026-06-24)
- **Python 3.14.6**

R packages: GSVA, DESeq2, GEOquery, org.Hs.eg.db, AnnotationDbi, car,
vegan, limma, ggplot2, dplyr, tidyr, patchwork, ggrepel, openxlsx,
pyreadr (via reticulate where used).

Python packages: numpy, pandas, scipy, statsmodels, matplotlib,
seaborn, pyreadr, gseapy.

Run `sessionInfo()` and `pip list` to record exact versions if
reproducing.

---

## Running the analysis

**Paths are hardcoded** to `/Users/mmani3/Downloads/Human Islets` in
most scripts. Change the `setwd()` line at the top of each before
running.

Suggested order:

```
discovery/00_inspect_download.R        structure of the downloaded files
discovery/01_map_and_detect.R          identifier mapping, detection frequency
discovery/02_model_crinophagy.R        pathway scoring, discovery models
discovery/03_pipeline.py               same analysis, independent Python implementation
discovery/04_robustness.py             permutation, bootstrap, leave-one-out
discovery/05_unc13d_check.py           UNC13D dropout diagnostics
discovery/06_rab3a_and_transcript.py   RAB3A robustness; transcript-level check
discovery/07_unc13_paralogues.py       Munc13 family, discovery proteomics

replication/08_gse164416_inspect.py    fetch and inspect GEO metadata
replication/09_replication_gse164416.R replication cohort analysis
replication/13_canonical_analysis.R    canonical dual-model analysis
replication/15_volcano_and_df_check.R  degrees-of-freedom sensitivity
replication/11_glycemic_continuum.R    IGT and T3cD comparisons

specificity/12_specificity_panel.R     pre-registered comparator panel
specificity/16_rab_layer_figure.R      Rab GTPase and GEF/GAP layer
specificity/17_exocytosis_panel.R      regulated exocytosis panel
specificity/18_munc13_figure.R         Munc13 paralogues
specificity/28_caps_family.R           CAPS family
specificity/25_extended_autophagy_panel.R   exploratory autophagy-lysosome panel
specificity/26_donor_acceptor_panel.R  donor-side vs acceptor-side panel

figures/...                            figure generation
tables/27_supplementary_tables.R       supplementary workbook
```

**`13_canonical_analysis.R` is the key script.** It builds `mat`,
`pd2`, `bg1` and `bg2`, the objects most later scripts depend on. The
session workspace (`session_gse164416.RData`) is not distributed and
must be rebuilt by running script 13, then saved:

```r
save(mat, pd2, bg1, bg2, ratio_all, file = "session_gse164416.RData")
```

---

## Notes on specific scripts

**`discovery/01_map_and_detect.R`** supersedes an earlier detection
script written before it was established that the two proteomics
matrices use different identifier systems (UniProt in the unprocessed
matrix, Entrez in the processed one). Gene sets are mapped separately
for each.

**Identifier handling.** All accessions are resolved programmatically
through `org.Hs.eg.db` rather than entered by hand. An incorrect
accession produces a silent absence rather than an error, and absence
is not distinguishable from biological non-detection without
verification.

**`replication/11_glycemic_continuum.R`** analyses IGT and T3cD donors.
These comparisons are reported in the manuscript but do not carry the
primary claim.

**`specificity/25_extended_autophagy_panel.R` and
`specificity/26_donor_acceptor_panel.R`** produced null results, which
are reported as such in the manuscript. The code is included because
the analyses were performed.

**`figures/22_transcript_protein_groups.R`** produced a panel that was
not used in the final manuscript: 19 of the 46 panel genes, including
STX7 and CTSD, were not quantified in both molecular layers, so the
figure could not display the principal candidates. Retained for
completeness.

**`discovery/05_unc13d_check.py`** addresses a detection confound in
the discovery proteomics. UNC13D was quantified in 97 of 129 donors
(75.2%), with detection differing by glycemic group. This motivated
the transcript-level replication, where no detection limit applies.

---

## Reproducibility

Random seeds are fixed for all permutation and bootstrap procedures.
Permutation P values may vary in the third decimal place between runs
and between the R and Python implementations, which use different
random number generators; the values reported in the manuscript derive
from a single frozen set of output files.

The primary variance analysis was implemented independently in R and
Python and gave identical results. This is reported as a check against
implementation error, not as independent evidence.

---

## Citation

If you use this code, please cite the manuscript (details to follow)
and the gene set deposit:

> Mani M, Martin TFJ. Pre-registered crinophagy and macroautophagy
> gene sets. Zenodo, 2026. https://doi.org/10.5281/zenodo.23060343

---

## License

[ MIT for code]

---

## Contact

Muralidharan Mani — mmani3@wisc.edu
