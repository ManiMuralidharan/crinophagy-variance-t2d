#!/usr/bin/env python3
"""
Robustness checks for the crinophagy variance-heterogeneity finding.

Observed (03_pipeline.py):
  crinophagy     Brown-Forsythe W = 4.144, p = 0.0182
  macroautophagy Brown-Forsythe W = 0.481, p = 0.6196
  SD ratio Type2 / NoDiabetes = 62.92 / 38.97 = 1.61
  n = 87 / 20 / 16 after covariate adjustment

The T2D arm has n = 16. Variance tests are less stable than mean tests
at that size, so the point estimate alone is not publishable. This
script attaches uncertainty and tests fragility.

  1. Permutation test          - distribution-free p for W
  2. Bootstrap CI              - uncertainty on the SD ratio
  3. Leave-one-out             - does one donor carry the result?
  4. Per-protein decomposition - is it the Ca2+ machinery specifically?

Usage:
  python 04_robustness.py
  python 04_robustness.py --nperm 50000 --seed 42
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

NUMERIC = ["crino_raw", "macro_raw", "crino_norm", "macro_norm",
           "crino_no_unc", "hba1c", "donorage", "bodymassindex",
           "predistributionculturetime", "percentieqrecoverypostculture",
           "pdinsulinperieq", "pdinsulindnaratio", "pdinsulincontent"]

NEED = ["crino_raw", "macro_raw", "donorage", "bodymassindex",
        "predistributionculturetime", "percentieqrecoverypostculture"]


def read_csv_safe(path):
    """pandas reads the string 'None' as NaN, which would silently drop
    all 91 non-diabetic donors. Never read this file any other way."""
    d = pd.read_csv(path, na_values=["NA", ""], keep_default_na=False)
    for c in NUMERIC:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d["diagnosis_computed"] = d["diagnosis_computed"].replace(
        "None", "NoDiabetes")
    return d


def adjusted_residuals(d):
    """Covariate-adjust once; all four checks operate on these."""
    sub = d.dropna(subset=NEED).copy()
    for y in ["crino_raw", "macro_raw"]:
        sub["r_" + y] = smf.ols(f"{y} ~ {COVARS}", data=sub).fit().resid
    return sub


def bf_stat(values, labels):
    groups = [values[labels == g] for g in GROUPS]
    groups = [g for g in groups if len(g) > 1]
    return stats.levene(*groups, center="median").statistic


# ---------------------------------------------------------------------
# 1. Permutation test
# ---------------------------------------------------------------------

def permutation_test(sub, col, nperm, rng):
    vals = sub["r_" + col].to_numpy()
    labs = sub["diagnosis_computed"].to_numpy()

    obs = bf_stat(vals, labs)
    null = np.empty(nperm)
    perm = labs.copy()
    for i in range(nperm):
        rng.shuffle(perm)
        null[i] = bf_stat(vals, perm)

    # +1 correction: a permutation p can never be exactly zero
    p = (np.sum(null >= obs) + 1) / (nperm + 1)
    return obs, p, null


# ---------------------------------------------------------------------
# 2. Bootstrap CI on the SD ratio
# ---------------------------------------------------------------------

def bootstrap_sd_ratio(sub, col, nboot, rng):
    a = sub.loc[sub.diagnosis_computed == "Type2", "r_" + col].to_numpy()
    b = sub.loc[sub.diagnosis_computed == "NoDiabetes", "r_" + col].to_numpy()

    obs = a.std(ddof=1) / b.std(ddof=1)

    boot = np.empty(nboot)
    for i in range(nboot):
        # resample within group, preserving group sizes
        aa = rng.choice(a, size=len(a), replace=True)
        bb = rng.choice(b, size=len(b), replace=True)
        boot[i] = aa.std(ddof=1) / bb.std(ddof=1)

    lo, hi = np.percentile(boot, [2.5, 97.5])
    return obs, lo, hi, boot


# ---------------------------------------------------------------------
# 3. Leave-one-out
# ---------------------------------------------------------------------

def leave_one_out(sub, col):
    """Drop each Type2 donor in turn. If p crosses 0.05 for any single
    exclusion, the finding rests on one sample and must be reported
    as such."""
    t2 = sub.index[sub.diagnosis_computed == "Type2"]
    rows = []
    for idx in t2:
        s = sub.drop(index=idx)
        groups = [s.loc[s.diagnosis_computed == g, "r_" + col] for g in GROUPS]
        lev = stats.levene(*groups, center="median")
        rows.append(dict(dropped=sub.loc[idx, "record_id"],
                         W=lev.statistic, p=lev.pvalue))
    return pd.DataFrame(rows).sort_values("p", ascending=False)


# ---------------------------------------------------------------------
# 4. Per-protein decomposition
# ---------------------------------------------------------------------

def per_protein(sub, expr_path, sets_path):
    """Where does the dispersion live? A diffuse pathway effect and a
    Ca2+-machinery-specific effect are different papers."""
    try:
        import pyreadr
    except ImportError:
        print("  pyreadr not installed; skipping. pip install pyreadr")
        return None

    prot = pyreadr.read_r(str(expr_path))[None]
    sets_df = pd.read_csv(sets_path, sep="\t")

    expr = prot.set_index(prot.columns[0]).astype(float)
    expr.index = expr.index.astype(str)

    keep = [r for r in sub.record_id if r in expr.columns]
    meta = sub.set_index("record_id").loc[keep]

    rows = []
    for _, g in sets_df.iterrows():
        eid = str(g.ENTREZID)
        if eid not in expr.index:
            continue
        v = expr.loc[eid, keep]
        df = meta.assign(y=v.to_numpy())
        df = df.dropna(subset=["y"])
        if df.y.std() == 0:
            continue

        # adjust each protein for the same covariates
        df["r"] = smf.ols(f"y ~ {COVARS}", data=df).fit().resid
        a = df.loc[df.diagnosis_computed == "Type2", "r"]
        b = df.loc[df.diagnosis_computed == "NoDiabetes", "r"]
        if len(a) < 3 or len(b) < 3 or b.std() == 0:
            continue

        # Levene on two groups, and the variance ratio with an F test
        lev = stats.levene(a, b, center="median")
        F = a.var(ddof=1) / b.var(ddof=1)
        rows.append(dict(symbol=g.SYMBOL, set=g.set,
                         sd_ratio=a.std(ddof=1) / b.std(ddof=1),
                         F=F, p=lev.pvalue))

    out = pd.DataFrame(rows)
    if len(out) == 0:
        return None

    # Benjamini-Hochberg within the crinophagy set
    from statsmodels.stats.multitest import multipletests
    out["q"] = np.nan
    for s in out.set.unique():
        m = out.set == s
        out.loc[m, "q"] = multipletests(out.loc[m, "p"], method="fdr_bh")[1]

    return out.sort_values("sd_ratio", ascending=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--nboot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--csv", default="modelling_dataset.csv")
    args = ap.parse_args()

    rng = np.random.default_rng(args.seed)

    d = read_csv_safe(BASE / args.csv)
    sub = adjusted_residuals(d)

    print(f"Donors after covariate adjustment: {len(sub)}")
    print(sub.diagnosis_computed.value_counts().reindex(GROUPS).to_string())

    # ---- 1 ----
    print(f"\n===== 1. PERMUTATION TEST ({args.nperm:,} shuffles) =====")
    for col in ["crino_raw", "macro_raw"]:
        obs, p, null = permutation_test(sub, col, args.nperm, rng)
        print(f"  {col:12s} W_obs={obs:6.3f}  perm p={p:.4f}  "
              f"(null mean W={null.mean():.3f})")

    # ---- 2 ----
    print(f"\n===== 2. BOOTSTRAP SD RATIO ({args.nboot:,} resamples) =====")
    for col in ["crino_raw", "macro_raw"]:
        obs, lo, hi, boot = bootstrap_sd_ratio(sub, col, args.nboot, rng)
        excl = "excludes 1.0" if lo > 1 or hi < 1 else "INCLUDES 1.0"
        print(f"  {col:12s} Type2/NoDiabetes SD ratio = {obs:.3f}  "
              f"95% CI [{lo:.3f}, {hi:.3f}]  {excl}")

    # ---- 3 ----
    print("\n===== 3. LEAVE-ONE-OUT (Type2 donors) =====")
    loo = leave_one_out(sub, "crino_raw")
    print(loo.round(4).to_string(index=False))
    worst = loo.iloc[0]
    print(f"\n  Worst case: dropping {worst.dropped} gives p = {worst.p:.4f}")
    print("  ROBUST" if worst.p < 0.05 else
          "  FRAGILE: one donor can overturn the result. Report this.")

    # ---- 4 ----
    print("\n===== 4. PER-PROTEIN DECOMPOSITION =====")
    pp = per_protein(sub,
                     BASE / "proteomics" / "proc_prot.rds",
                     BASE / "frozen_gene_sets_v2_clean.txt")
    if pp is not None:
        print("\n  Top 12 by SD ratio (Type2 / NoDiabetes):")
        print(pp.head(12).round(4).to_string(index=False))

        core = pp[pp.symbol.isin(["UNC13D", "CADPS", "RAB3A", "RAB3GAP1"])]
        print("\n  Ca2+ machinery:")
        print(core.round(4).to_string(index=False))

        print("\n  Median SD ratio by set:")
        print(pp.groupby("set").sd_ratio.median().round(3).to_string())

        pp.to_csv(BASE / "per_protein_variance.csv", index=False)
        print("\n  Wrote per_protein_variance.csv")


if __name__ == "__main__":
    main()
