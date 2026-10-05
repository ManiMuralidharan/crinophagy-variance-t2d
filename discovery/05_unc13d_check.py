#!/usr/bin/env python3
"""
UNC13D dropout diagnostics.

Context (04_robustness.py, per-protein decomposition):
  UNC13D  SD ratio 2.07  F=4.27  Levene p=0.0017  q=0.022
  RAB3A   SD ratio 1.73  F=2.98  Levene p=0.0014  q=0.022
  (only two proteins clearing FDR across the 31-protein set)

THE PROBLEM: UNC13D was quantified in 97/129 donors (75.2%), while
almost every other set member sat at 100%. Elevated variance in a
protein with more missing data is exactly what technical dropout
produces. Until that is excluded, the UNC13D result is not safe.
RAB3A at 100% detection has no such vulnerability.

Tests run here:
  1. Is UNC13D missingness differential by glycemic group?
  2. Does the variance ratio survive restriction to complete cases?
  3. Does UNC13D residual variance track MS depth per donor?
  4. Is low UNC13D confounded with overall low signal in that donor?
  5. Does the pathway test sharpen when reduced to UNC13D + RAB3A?

Usage:
  python 05_unc13d_check.py
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

BASE = Path("/Users/mmani3/Downloads/Human Islets")

COVARS = ("donorage + C(donorsex) + bodymassindex "
          "+ predistributionculturetime + percentieqrecoverypostculture")

GROUPS = ["NoDiabetes", "Pre.T2D", "Type2"]

NUMERIC = ["crino_raw", "macro_raw", "crino_no_unc", "hba1c", "donorage",
           "bodymassindex", "predistributionculturetime",
           "percentieqrecoverypostculture", "pdinsulinperieq",
           "pdinsulindnaratio", "pdinsulincontent"]

NEED = ["crino_raw", "donorage", "bodymassindex",
        "predistributionculturetime", "percentieqrecoverypostculture"]

UNC13D, RAB3A = "201294", "5864"      # Entrez IDs


def read_csv_safe(path):
    """pandas reads 'None' as NaN, which would drop all 91 non-diabetic
    donors. Never read this file any other way."""
    d = pd.read_csv(path, na_values=["NA", ""], keep_default_na=False)
    for c in NUMERIC:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d["diagnosis_computed"] = d["diagnosis_computed"].replace(
        "None", "NoDiabetes")
    return d


def load_matrices():
    import pyreadr
    proc = pyreadr.read_r(str(BASE / "proteomics" / "proc_prot.rds"))[None]
    raw = pyreadr.read_r(str(BASE / "proteomics" / "unproc_prot.rds"))[None]

    proc = proc.set_index(proc.columns[0])
    proc.index = proc.index.astype(str)
    proc = proc.astype(float)

    raw = raw.set_index(raw.columns[0])
    raw.index = raw.index.astype(str)
    raw = raw.astype(float)
    return proc, raw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="modelling_dataset.csv")
    args = ap.parse_args()

    d = read_csv_safe(BASE / args.csv)
    proc, raw = load_matrices()

    # UNC13D UniProt accession, since the raw matrix is UniProt-keyed
    unc_uni = "Q70J99"

    # ---- 1. Differential missingness ---------------------------------
    print("===== 1. IS UNC13D MISSINGNESS DIFFERENTIAL BY GROUP? =====")

    if unc_uni in raw.index:
        v = raw.loc[unc_uni]
        det = pd.DataFrame({
            "record_id": v.index,
            "detected": (~v.isna() & (v != 0)).to_numpy()
        })
        m = det.merge(d[["record_id", "diagnosis_computed"]], on="record_id")

        tab = pd.crosstab(m.diagnosis_computed, m.detected)
        tab = tab.reindex(GROUPS)
        print(tab.to_string())
        print("\n  Detection rate by group:")
        for g in GROUPS:
            r = m[m.diagnosis_computed == g]
            print(f"    {g:12s} {r.detected.sum():3d}/{len(r):3d} "
                  f"= {100*r.detected.mean():5.1f}%")

        chi = stats.chi2_contingency(tab.to_numpy())
        print(f"\n  Chi-square p = {chi[1]:.4f}")
        print("  If significant, dropout is confounded with disease status")
        print("  and the variance result cannot be interpreted as biology.")
    else:
        print(f"  {unc_uni} not in raw matrix; check the accession.")

    # ---- 2. Complete-case variance ratio ------------------------------
    print("\n===== 2. VARIANCE RATIO IN COMPLETE CASES ONLY =====")

    sub = d.dropna(subset=NEED).copy()
    keep = [r for r in sub.record_id if r in proc.columns]
    meta = sub.set_index("record_id").loc[keep]

    for name, eid in [("UNC13D", UNC13D), ("RAB3A", RAB3A)]:
        if eid not in proc.index:
            print(f"  {name}: {eid} not in processed matrix")
            continue

        df = meta.assign(y=proc.loc[eid, keep].to_numpy()).dropna(subset=["y"])
        df["r"] = smf.ols(f"y ~ {COVARS}", data=df).fit().resid

        a = df.loc[df.diagnosis_computed == "Type2", "r"]
        b = df.loc[df.diagnosis_computed == "NoDiabetes", "r"]
        lev = stats.levene(a, b, center="median")

        # bootstrap the ratio
        rng = np.random.default_rng(1)
        boot = np.array([
            rng.choice(a, len(a), True).std(ddof=1) /
            rng.choice(b, len(b), True).std(ddof=1)
            for _ in range(5000)])
        lo, hi = np.percentile(boot, [2.5, 97.5])

        print(f"  {name:7s} n(T2D)={len(a):3d} n(ND)={len(b):3d}  "
              f"SD ratio={a.std(ddof=1)/b.std(ddof=1):.3f}  "
              f"95% CI [{lo:.3f}, {hi:.3f}]  Levene p={lev.pvalue:.4f}")

    # ---- 3. Does variance track MS depth? -----------------------------
    print("\n===== 3. UNC13D RESIDUAL vs MS DEPTH =====")

    # depth proxy: proteins quantified per donor in the RAW matrix
    depth = (~raw.isna() & (raw != 0)).sum(axis=0)
    depth.name = "n_proteins"
    dep = depth.reset_index()
    dep.columns = ["record_id", "n_proteins"]

    df = meta.reset_index().merge(dep, on="record_id")
    df["y"] = proc.loc[UNC13D, df.record_id].to_numpy()
    df = df.dropna(subset=["y", "n_proteins"])

    print(f"  Proteins quantified per donor: "
          f"median={df.n_proteins.median():.0f}, "
          f"range {df.n_proteins.min():.0f}-{df.n_proteins.max():.0f}")

    r, p = stats.pearsonr(df.n_proteins, df.y)
    print(f"  UNC13D abundance vs depth: r={r:.3f}, p={p:.4f}")

    df["r_unc"] = smf.ols(f"y ~ {COVARS}", data=df).fit().resid
    r2, p2 = stats.pearsonr(df.n_proteins, df.r_unc.abs())
    print(f"  |residual| vs depth:       r={r2:.3f}, p={p2:.4f}")
    print("  A strong negative correlation would mean the extra variance")
    print("  is measurement noise in shallow runs, not biology.")

    # depth as an extra covariate
    m = smf.ols(f"y ~ n_proteins + {COVARS}", data=df).fit()
    df["r_adj"] = m.resid
    a = df.loc[df.diagnosis_computed == "Type2", "r_adj"]
    b = df.loc[df.diagnosis_computed == "NoDiabetes", "r_adj"]
    lev = stats.levene(a, b, center="median")
    print(f"\n  After adjusting for depth: SD ratio="
          f"{a.std(ddof=1)/b.std(ddof=1):.3f}, Levene p={lev.pvalue:.4f}")

    # ---- 4. Is depth itself different by group? -----------------------
    print("\n===== 4. IS MS DEPTH DIFFERENT BY GROUP? =====")
    for g in GROUPS:
        s = df.loc[df.diagnosis_computed == g, "n_proteins"]
        print(f"  {g:12s} n={len(s):3d}  median={s.median():.0f}  "
              f"SD={s.std():.0f}")
    kw = stats.kruskal(*[df.loc[df.diagnosis_computed == g, "n_proteins"]
                         for g in GROUPS])
    print(f"  Kruskal-Wallis p = {kw.pvalue:.4f}")

    # ---- 5. Reduced two-protein signature -----------------------------
    print("\n===== 5. REDUCED SIGNATURE: UNC13D + RAB3A =====")

    z = proc.loc[[UNC13D, RAB3A], keep]
    z = z.sub(z.mean(axis=1), axis=0).div(z.std(axis=1), axis=0)
    df2 = meta.assign(y=z.mean(axis=0).to_numpy()).dropna(subset=["y"])
    df2["r"] = smf.ols(f"y ~ {COVARS}", data=df2).fit().resid

    groups = [df2.loc[df2.diagnosis_computed == g, "r"] for g in GROUPS]
    lev = stats.levene(*groups, center="median")
    print(f"  Brown-Forsythe (3 groups): W={lev.statistic:.3f}  "
          f"p={lev.pvalue:.4f}")
    for g, s in zip(GROUPS, groups):
        print(f"    {g:12s} n={len(s):3d}  SD={s.std():.3f}")

    # leave-one-out on the reduced signature
    t2 = df2.index[df2.diagnosis_computed == "Type2"]
    ps = []
    for idx in t2:
        s = df2.drop(index=idx)
        gg = [s.loc[s.diagnosis_computed == g, "r"] for g in GROUPS]
        ps.append(stats.levene(*gg, center="median").pvalue)
    print(f"\n  Leave-one-out worst p = {max(ps):.4f}")
    print("  ROBUST" if max(ps) < 0.05 else
          "  Still fragile, but compare against the 31-protein version.")


if __name__ == "__main__":
    main()
