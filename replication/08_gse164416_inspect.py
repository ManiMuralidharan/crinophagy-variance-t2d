#!/usr/bin/env python3
"""
Replication cohort: GSE164416 (Wigger, Barovic & Brunner et al.,
Nat Metab 2021). LCM-isolated islets from pancreatectomised living
donors across the glycemic continuum: ND, IGT, T3cD, T2D.

WHY THIS COHORT. The paper's own headline is that diabetic islets are
heterogeneous where non-diabetic ones are not, so elevated
inter-individual variance in T2D has independent precedent. The
question here is narrow: does RAB3A specifically show it?

WHAT DIFFERS FROM HumanIslets. Surgical LCM rather than enzymatically
isolated cultured islets; RNA-seq rather than MS; different donors;
no culture time to adjust for. Replication across that much
methodological difference would be meaningful. Note that transcript
variance and protein variance are different quantities: a null here
does not refute the protein finding, it just fails to support it.

STAGE 1 (this script): download and inspect. Do not analyse yet.
The covariate set and group labels depend on what the metadata
actually contains, which we cannot assume.

Usage:
  cd "/Users/mmani3/Downloads/Human Islets"
  python3 08_gse164416_inspect.py
"""

import gzip
import io
import re
import tarfile
import urllib.request
from pathlib import Path

import pandas as pd

BASE = Path("/Users/mmani3/Downloads/Human Islets")
REP = BASE / "replication_GSE164416"
REP.mkdir(exist_ok=True)

GEO = "GSE164416"
SERIES_MATRIX = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/"
    f"{GEO}/matrix/{GEO}_series_matrix.txt.gz"
)
SUPP_DIR = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE164nnn/"
    f"{GEO}/suppl/"
)


def fetch(url, dest):
    if dest.exists():
        print(f"  already have {dest.name}")
        return dest
    print(f"  downloading {url}")
    urllib.request.urlretrieve(url, dest)
    print(f"  -> {dest.name} ({dest.stat().st_size/1e6:.1f} MB)")
    return dest


def main():
    print("=" * 64)
    print("STAGE 1: fetch and inspect GSE164416")
    print("=" * 64)

    # ---- 1. Series matrix: sample metadata -------------------------
    print("\n--- Series matrix (metadata) ---")
    sm = REP / f"{GEO}_series_matrix.txt.gz"
    try:
        fetch(SERIES_MATRIX, sm)
    except Exception as e:
        print(f"  FAILED: {e}")
        print("  Fetch manually from:")
        print(f"  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GEO}")
        return

    # The series matrix header holds !Sample_* lines. Characteristics
    # lines carry the clinical variables we need.
    meta_lines = []
    with gzip.open(sm, "rt", errors="replace") as f:
        for line in f:
            if line.startswith("!series_matrix_table_begin"):
                break
            meta_lines.append(line.rstrip("\n"))

    print(f"\n  {len(meta_lines)} header lines")
    print("\n  Sample-level fields present:")
    for ln in meta_lines:
        if ln.startswith("!Sample_"):
            key = ln.split("\t")[0]
            n_vals = len(ln.split("\t")) - 1
            preview = ln.split("\t")[1][:60] if n_vals else ""
            print(f"    {key:38s} n={n_vals:4d}  {preview}")

    # Parse characteristics into a tidy frame
    rows = {}
    for ln in meta_lines:
        parts = ln.split("\t")
        key = parts[0].strip('!"')
        vals = [p.strip('"') for p in parts[1:]]
        if key.startswith("Sample_characteristics") or \
           key in ("Sample_title", "Sample_geo_accession",
                   "Sample_source_name_ch1"):
            # characteristics come as "label: value"
            if vals and ":" in vals[0]:
                label = vals[0].split(":")[0].strip()
                rows[label] = [v.split(":", 1)[1].strip()
                               if ":" in v else v for v in vals]
            else:
                rows[key] = vals

    md = pd.DataFrame(rows)
    md.to_csv(REP / "gse164416_metadata.csv", index=False)

    print(f"\n  Parsed metadata: {md.shape}")
    print(f"  Columns: {list(md.columns)}")
    print("\n  First 3 samples:")
    print(md.head(3).to_string())

    # What are the group labels actually called?
    print("\n  Value counts for likely grouping columns:")
    for c in md.columns:
        u = md[c].nunique()
        if 2 <= u <= 8:
            print(f"\n    {c}:")
            print("      " + md[c].value_counts().to_string()
                  .replace("\n", "\n      "))

    # ---- 2. Supplementary files: the count matrix ------------------
    print("\n--- Supplementary files (expression) ---")
    print("  GEO supplementary listing:")
    print(f"  {SUPP_DIR}")
    try:
        with urllib.request.urlopen(SUPP_DIR) as r:
            html = r.read().decode("utf-8", errors="replace")
        files = sorted(set(re.findall(r'href="([^"]+\.(?:txt|csv|tsv|xlsx|gz|tar))"',
                                      html)))
        for f in files:
            print(f"    {f}")
        print("\n  Download the counts file with:")
        print(f"    curl -O {SUPP_DIR}<filename>")
    except Exception as e:
        print(f"  Could not list: {e}")
        print(f"  Browse manually:")
        print(f"  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={GEO}")

    # ---- 3. What we need next --------------------------------------
    print("\n" + "=" * 64)
    print("WHAT TO REPORT BACK")
    print("=" * 64)
    print("""
  1. The grouping column name and its levels (ND / IGT / T3cD / T2D
     or whatever GEO actually calls them), with n per group.

  2. Which covariates exist: age, sex, BMI, HbA1c. There is no
     culture time here, which simplifies the model, but check what
     surgical or batch variables are recorded instead.

  3. The supplementary filename holding the count matrix, and
     whether rows are Ensembl IDs, gene symbols, or Entrez.

  4. The T2D group size. This is the number that decides whether a
     variance test is worth running at all. Below ~12 it will not
     be informative, and we should say so rather than run it.

  Then I can write the analysis script against the real structure,
  as we did with the HumanIslets download.
""")

    print(f"  Metadata saved: {REP / 'gse164416_metadata.csv'}")


if __name__ == "__main__":
    main()
