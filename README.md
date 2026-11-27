# crinophagy-variance-t2d

Code, frozen gene sets and summary results for:

> Mani M., et al. *Inter-individual variability of crinophagy-related molecules in human islets from donors with type 2 diabetes.* International Journal of Molecular Sciences (manuscript IJMS-4686697, revised version).

The study is exploratory. It examines inter-donor dispersion (residual SD ratios) of a frozen 31-member crinophagy-related gene set and a 14-member macroautophagy comparator in HumanIslets discovery proteomics and an independent bulk RNA-seq cohort (GSE164416).

## Two analysis stages

This repository contains two clearly separated stages. Please keep them distinct when citing numbers.

| Stage | Folder | When | What it is |
|---|---|---|---|
| **Historical analysis** | `scripts/historical/`, `results/historical/` | September 2026 | The analyses reported in the originally submitted manuscript. Kept unchanged for transparency. |
| **Revision analyses** | `scripts/revision/`, `results/revision/` | November 2026 | Checks, corrections and sensitivity analyses added in response to peer review. |

Some historical statements were corrected during revision. The main corrections are:

- Historical per-protein SD ratios were computed from the **processed** proteomics matrix (26 of 31 crinophagy genes), not the unprocessed matrix as originally stated. A separate raw-matrix sensitivity (all 31 genes) is in `results/revision/41_*` and `45_*`.
- The gene sets were frozen (`genesets/`), but the statistical models, candidate hierarchy and testing families were **not** preregistered. All analyses are exploratory.
- The Figure 1C three-group permutation P value is **0.013** (recomputed from `04_robustness.py`), not 0.017.
- Discovery donor deletion: **three of 16** T2D donors raise the omnibus P above 0.05 (fixed residuals; `figures/` and `results/historical/C_dispersion_loo.csv`).
- Multiplicity: BH correction is reported within the 31-gene and 14-gene families, jointly across 45 genes, and Holm across UNC13D/RAB3A (`results/revision/01_replication_multiplicity.csv`).
- The extended exploratory panel contains **93** genuinely additional genes, not 107 (`results/historical/C16_*`).
- Figure 2C was regenerated with 1,000,000 random gene sets and Figure S4 without a filter that dropped negative processed means (`figures/revision/`).

## Repository layout

```
genesets/                 frozen gene sets v1, v2, v2_clean, v3 (final 31 + 14) and identifier mapping
scripts/historical/       00-36: original discovery, replication, panel and figure scripts
scripts/revision/         40-77: revision checks and sensitivity analyses
results/historical/       summary tables from the original analysis (C_*.csv, replication tables)
results/revision/         summary tables and audit notes from the revision
figures/                  original figure PDFs
figures/revision/         regenerated Figure 2C and Figure S4 (PDF and editable SVG)
session_info/             R and Python environment records
```

## Data access

Raw and donor-level data are **not** included.

- **Discovery cohort:** HumanIslets.com (registration required). Files used: `proteomics/proc_prot.rds`, `proteomics/unproc_prot.rds`, `Bulk RNA/proc_rnaseq.rds`, `donor.csv`.
- **Independent cohort:** GEO accession [GSE164416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE164416) (HTSeq counts and series metadata).
- **Single-cell context:** published supplementary tables of Lawlor et al. (GEO GSE86473); no reanalysis.

## How to run

1. Download the data above into one folder, keeping the subfolders `proteomics/`, `Bulk RNA/` and `replication_GSE164416/`.
2. Point the scripts at that folder:
   ```bash
   export HUMAN_ISLETS_DIR=/path/to/data/folder
   ```
3. Run scripts in numerical order within each stage. Python scripts take the data folder and an output folder as arguments, for example:
   ```bash
   python3 scripts/revision/74_zscore_pathway_score.py "$HUMAN_ISLETS_DIR" out/
   ```

Software: R 4.6.1 (GSVA, DESeq2 1.52.0, vegan, car 3.1-5, org.Hs.eg.db 3.23.1) and Python 3 (NumPy, pandas, SciPy, statsmodels, pyreadr). Exact versions are in `session_info/`; historical analyses used GSVA 2.6.6 and vegan 2.7.5.

## Script guide

**Historical (00-36)**

| Scripts | Purpose |
|---|---|
| 00-02 | Data inspection, gene mapping, discovery pathway scores and models |
| 03-07 | Python pipeline, robustness (permutation and bootstrap), candidate checks |
| 08-13 | GSE164416 replication, high-variance enrichment, Model 1/Model 2 analysis |
| 14-30 | Figures, panels (exocytosis, Munc13, CAPS, extended autophagy-lysosome), supplementary tables |
| 31-36 | Table rebuilding, INS sensitivity, reconciliation, Figure 5A template |

Notes: `31_rebuild_supplementary_tables_fixed.R` supersedes the unfixed version. `36_fig5A.R` is a template whose discovery intervals are supplied by `results/revision/75_discovery_fig5A_intervals.csv`.

**Revision (40-77)**

| Scripts | Purpose | Key outputs |
|---|---|---|
| 40-45 | Disease-label permutation, joint bootstrap, resampling audit, multiplicity families | `C1_C2_between_set.csv`, `C9_multiplicity_both_panels.csv` |
| 50-57 | Discovery checks, cross-cohort comparison, dispersion models, calibration, pathway CI | `C_nested_dispersion.csv`, `C_calibration.csv`, `C_pathway_ci_vs_bf.csv` |
| 59-60 | Post-hoc fusion-subset sensitivity and independent-cohort validation | `results/revision/01-07` |
| 61-66 | Paired protein/RNA analysis, coverage, interval coverage | `19_same_donor_covariate_adjusted_SDR.csv`, `15_*`, `16_*` |
| 65, 67-72 | Raw-matrix audit, UNC13D missingness, preprocessing provenance, five excluded genes, raw 31-gene sensitivity | `38_*` to `47_*` |
| 74 | Standardized-abundance (mean z-score) pathway score, both cohorts | `74_zscore_pathway_score_results.csv` |
| 75 | Recomputed Figure 5A discovery intervals | `75_discovery_fig5A_intervals.csv` |
| 76/77 | Regenerated Figure 2C (1,000,000 random sets) and Figure S4 | `76_*.csv`, `77_*.csv`, `figures/revision/` |

Superseded versions are not included: `43_multiplicity_fixed.R`, `65_raw_proteomics_uniprot_audit_fixed.R` and `70_discovery_five_excluded_crinophagy_genes_FIXED.R` replace their unfixed versions.

## Citation and archive

Frozen gene sets: https://doi.org/10.5281/zenodo.23060343
Code and results (this release, v1.0-revision): https://doi.org/10.5281/zenodo.23264616
