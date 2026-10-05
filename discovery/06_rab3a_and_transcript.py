#!/usr/bin/env python3
"""
Two analyses, following from 05_unc13d_check.py.

PART A - RAB3A alone.
  RAB3A is detected in 129/129 donors, so it faces none of the dropout
  objections that compromised UNC13D. Observed:
    SD ratio (Type2/NoDiabetes) = 1.725, 95% CI [1.092, 2.276]
    Levene p = 0.0014
  This runs the full robustness suite on the single protein:
  permutation, bootstrap, leave-one-out, and the macroautophagy
  comparator proteins for context.

PART B - the UNC13D detection puzzle.
  UNC13D detection: 68.1% in non-diabetic, 88.9% T2D, 95.0% pre-T2D
  (chi-square p = 0.0146). That confounds its variance result. But it
  raises a real question: is Munc13-4 genuinely higher in diabetic
  islets, crossing the MS detection limit more often? The transcript
  layer has no dropout, so it can answer this cleanly in the 87 donors
  with both proteomics and bulk RNA-seq.

Usage:
  python 06_rab3a_and_transcript.py
  python 06_rab3a_and_transcript.py --nperm 50000
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

RAB3A, UNC13D = "5864", "201294"          # Entrez
RAB3A_UNI, UNC13D_UNI = "P20336", "Q70J99"

# macroautophagy proteins, as single-protein negative controls
MACRO = {"ATG3": "64422", "ATG5": "9474", "ATG7": "10533",
         "BECN1": "8678", "ULK1": "8408", "WIPI2": "26100",
         "MAP1LC3A": "84557", "GABARAP": "11337",
         "GABARAPL1": "23710", "GABARAPL2": "11345",
         "SQSTM1": "8878", "NBR1": "4077", "CALCOCO2": "10241"}


def read_csv_safe(path):
    """pandas reads 'None' as NaN, dropping all 91 non-diabetic donors."""
    d = pd.read_csv(path, na_values=["NA", ""], keep_default_na=False)
    for c in NUMERIC:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d["diagnosis_computed"] = d["diagnosis_computed"].replace(
        "None", "NoDiabetes")
    return d


def residualise(meta, values, keep):
    df = meta.assign(y=values).dropna(subset=["y"])
    df["r"] = smf.ols(f"y ~ {COVARS}", data=df).fit().resid
    return df


def bf_stat(vals, labs):
    g = [vals[labs == k] for k in GROUPS]
    g = [x for x in g if len(x) > 1]
    return stats.levene(*g, center="median").statistic


def part_a(d, proc, nperm, nboot, rng):
    print("=" * 62)
    print("PART A: RAB3A alone")
    print("=" * 62)

    sub = d.dropna(subset=NEED).copy()
    keep = [r for r in sub.record_id if r in proc.columns]
    meta = sub.set_index("record_id").loc[keep]

    df = residualise(meta, proc.loc[RAB3A, keep].to_numpy(), keep)

    print(f"\nDonors: {len(df)}")
    for g in GROUPS:
        s = df.loc[df.diagnosis_computed == g, "r"]
        print(f"  {g:12s} n={len(s):3d}  SD={s.std(ddof=1):10.2f}")

    # --- Brown-Forsythe, 3 groups ---
    groups = [df.loc[df.diagnosis_computed == g, "r"] for g in GROUPS]
    lev = stats.levene(*groups, center="median")
    print(f"\n  Brown-Forsythe: W={lev.statistic:.3f}  p={lev.pvalue:.4f}")

    # --- permutation ---
    vals = df["r"].to_numpy()
    labs = df["diagnosis_computed"].to_numpy()
    obs = bf_stat(vals, labs)
    null = np.empty(nperm)
    perm = labs.copy()
    for i in range(nperm):
        rng.shuffle(perm)
        null[i] = bf_stat(vals, perm)
    pperm = (np.sum(null >= obs) + 1) / (nperm + 1)
    print(f"  Permutation ({nperm:,}): p={pperm:.4f}")

    # --- bootstrap SD ratio ---
    a = df.loc[df.diagnosis_computed == "Type2", "r"].to_numpy()
    b = df.loc[df.diagnosis_computed == "NoDiabetes", "r"].to_numpy()
    boot = np.array([
        rng.choice(a, len(a), True).std(ddof=1) /
        rng.choice(b, len(b), True).std(ddof=1) for _ in range(nboot)])
    lo, hi = np.percentile(boot, [2.5, 97.5])
    ratio = a.std(ddof=1) / b.std(ddof=1)
    print(f"  SD ratio T2D/ND = {ratio:.3f}  95% CI [{lo:.3f}, {hi:.3f}]  "
          f"{'excludes 1.0' if lo > 1 else 'INCLUDES 1.0'}")

    # --- leave-one-out ---
    print("\n  Leave-one-out (Type2 donors):")
    t2 = df.index[df.diagnosis_computed == "Type2"]
    rows = []
    for idx in t2:
        s = df.drop(index=idx)
        gg = [s.loc[s.diagnosis_computed == g, "r"] for g in GROUPS]
        rows.append(dict(dropped=idx,
                         p=stats.levene(*gg, center="median").pvalue))
    loo = pd.DataFrame(rows).sort_values("p", ascending=False)
    print(loo.head(5).round(4).to_string(index=False))
    worst = loo.p.max()
    print(f"\n  Worst-case p = {worst:.4f}  ->  "
          f"{'ROBUST' if worst < 0.05 else 'FRAGILE'}")
    print(f"  (31-protein score worst case was 0.104; "
          f"UNC13D+RAB3A composite was 0.156)")

    # --- negative controls: same test on each macroautophagy protein ---
    print("\n  Negative controls (macroautophagy, single proteins):")
    ctrl = []
    for sym, eid in MACRO.items():
        if eid not in proc.index:
            continue
        c = residualise(meta, proc.loc[eid, keep].to_numpy(), keep)
        aa = c.loc[c.diagnosis_computed == "Type2", "r"]
        bb = c.loc[c.diagnosis_computed == "NoDiabetes", "r"]
        if len(aa) < 3 or bb.std() == 0:
            continue
        ctrl.append(dict(symbol=sym,
                         sd_ratio=aa.std(ddof=1) / bb.std(ddof=1),
                         p=stats.levene(aa, bb, center="median").pvalue))
    ctrl = pd.DataFrame(ctrl).sort_values("sd_ratio", ascending=False)
    print(ctrl.round(4).to_string(index=False))
    n_sig = (ctrl.p < 0.05).sum()
    print(f"\n  {n_sig}/{len(ctrl)} controls reach p<0.05 "
          f"(expect ~{0.05*len(ctrl):.1f} by chance)")
    print(f"  RAB3A sd_ratio {ratio:.3f} vs control max "
          f"{ctrl.sd_ratio.max():.3f}")

    return df


def part_b(d, proc, raw):
    print("\n" + "=" * 62)
    print("PART B: is UNC13D genuinely higher in diabetic islets?")
    print("=" * 62)

    import pyreadr
    rna = pyreadr.read_r(str(BASE / "Bulk RNA" / "proc_rnaseq.rds"))[None]

    # proc_rnaseq has 5 metadata columns: accession, gene_id, symbol,
    # name, biotype. Locate UNC13D by symbol, not by position.
    meta_cols = list(rna.columns[:5])
    print(f"\n  RNA metadata columns: {meta_cols}")

    sym_col = next((c for c in meta_cols
                    if rna[c].astype(str).str.contains("UNC13D").any()), None)
    if sym_col is None:
        print("  UNC13D not found in RNA metadata; skipping Part B.")
        return

    row = rna[rna[sym_col].astype(str) == "UNC13D"]
    if len(row) == 0:
        print("  UNC13D absent from the RNA matrix.")
        return

    expr = row.iloc[0, 5:].astype(float)
    expr.index = expr.index.astype(str)

    md = d.set_index("record_id")
    shared = [r for r in expr.index if r in md.index]
    print(f"  Donors with RNA and metadata: {len(shared)}")

    df = md.loc[shared].copy()
    df["unc_rna"] = expr.loc[shared].to_numpy()
    df = df.dropna(subset=["unc_rna"] + NEED[1:])

    print("\n  UNC13D mRNA by group:")
    for g in GROUPS:
        s = df.loc[df.diagnosis_computed == g, "unc_rna"]
        if len(s):
            print(f"    {g:12s} n={len(s):3d}  median={s.median():10.3f}  "
                  f"SD={s.std():8.3f}")

    # mean difference, adjusted
    m = smf.ols(f"unc_rna ~ C(diagnosis_computed) + {COVARS}",
                data=df).fit()
    for t in m.params.index:
        if "diagnosis" in t:
            lab = t.split("[T.")[-1].rstrip("]")
            print(f"    {lab:12s} beta={m.params[t]:8.3f}  "
                  f"p={m.pvalues[t]:.4f}")

    # continuous HbA1c
    m2 = smf.ols(f"unc_rna ~ hba1c + {COVARS}", data=df).fit()
    print(f"\n  vs HbA1c: beta={m2.params['hba1c']:.4f}, "
          f"p={m2.pvalues['hba1c']:.4f}, n={int(m2.nobs)}")

    # variance at transcript level, where there is no dropout
    df["r"] = smf.ols(f"unc_rna ~ {COVARS}", data=df).fit().resid
    groups = [df.loc[df.diagnosis_computed == g, "r"] for g in GROUPS]
    groups = [g for g in groups if len(g) > 1]
    lev = stats.levene(*groups, center="median")
    print(f"  Variance (Brown-Forsythe): W={lev.statistic:.3f}  "
          f"p={lev.pvalue:.4f}")

    # does protein detection predict transcript level?
    if UNC13D_UNI in raw.index:
        v = raw.loc[UNC13D_UNI]
        det = pd.Series((~v.isna() & (v != 0)).to_numpy(),
                        index=v.index.astype(str))
        common = [r for r in df.index if r in det.index]
        dd = df.loc[common].assign(detected=det.loc[common].to_numpy())
        hi = dd.loc[dd.detected, "unc_rna"]
        lo = dd.loc[~dd.detected, "unc_rna"]
        if len(lo) > 1:
            u = stats.mannwhitneyu(hi, lo)
            print(f"\n  mRNA where protein detected  (n={len(hi)}): "
                  f"median={hi.median():.3f}")
            print(f"  mRNA where protein missing   (n={len(lo)}): "
                  f"median={lo.median():.3f}")
            print(f"  Mann-Whitney p = {u.pvalue:.4f}")
            print("\n  If detected donors have higher mRNA, the dropout is")
            print("  abundance-driven, i.e. a detection-limit effect rather")
            print("  than random technical loss. That makes the group")
            print("  difference in detection itself interpretable.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nperm", type=int, default=10000)
    ap.add_argument("--nboot", type=int, default=10000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--csv", default="modelling_dataset.csv")
    args = ap.parse_args()

    import pyreadr
    rng = np.random.default_rng(args.seed)

    d = read_csv_safe(BASE / args.csv)

    proc = pyreadr.read_r(str(BASE / "proteomics" / "proc_prot.rds"))[None]
    proc = proc.set_index(proc.columns[0])
    proc.index = proc.index.astype(str)
    proc = proc.astype(float)

    raw = pyreadr.read_r(str(BASE / "proteomics" / "unproc_prot.rds"))[None]
    raw = raw.set_index(raw.columns[0])
    raw.index = raw.index.astype(str)
    raw = raw.astype(float)

    part_a(d, proc, args.nperm, args.nboot, rng)
    part_b(d, proc, raw)


if __name__ == "__main__":
    main()
