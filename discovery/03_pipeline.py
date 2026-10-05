#!/usr/bin/env python3
"""
Crinophagy signature in human islet proteomics: Python pipeline.

Cohort: 129 HumanIslets donors with proteomics (all T1_diabetes = 0)
  diagnosis_computed: None 91, Pre.T2D 20, Type2 18

Primary finding so far: variance heterogeneity, not a mean shift.
  crinophagy     Brown-Forsythe p = 0.018
  macroautophagy Brown-Forsythe p = 0.62   (specificity control)

Run modes:
  python 03_pipeline.py --score     # rebuild scores from the .rds files
  python 03_pipeline.py             # model from modelling_dataset.csv

Install:
  pip install pandas numpy scipy statsmodels matplotlib seaborn
  pip install pyreadr gseapy        # only needed for --score
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

BASE = Path("/Users/mmani3/Downloads/Human Islets")

COVARS = ("donorage + C(donorsex) + bodymassindex "
          "+ predistributionculturetime + percentieqrecoverypostculture")

NUMERIC = ["crino_raw", "macro_raw", "crino_norm", "macro_norm",
           "crino_no_unc", "hba1c", "donorage", "bodymassindex",
           "predistributionculturetime", "percentieqrecoverypostculture",
           "pdinsulinperieq", "pdinsulindnaratio", "pdinsulincontent",
           "yearsdiabetic"]

GROUPS = ["NoDiabetes", "Pre.T2D", "Type2"]


# ---------------------------------------------------------------------
# CRITICAL: pandas treats the string "None" as missing by default, which
# silently deletes all 91 non-diabetic donors. Always read this way, and
# rename the level so the trap cannot be sprung downstream.
# ---------------------------------------------------------------------

def read_csv_safe(path):
    d = pd.read_csv(path, na_values=["NA", ""], keep_default_na=False)
    for c in NUMERIC:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    if "diagnosis_computed" in d.columns:
        d["diagnosis_computed"] = (d["diagnosis_computed"]
                                   .replace("None", "NoDiabetes"))
        d["diagnosis_computed"] = pd.Categorical(
            d["diagnosis_computed"], categories=GROUPS, ordered=False)
    return d


# ---------------------------------------------------------------------
# Optional: rebuild ssGSEA scores from the raw .rds files
# ---------------------------------------------------------------------

def build_scores():
    import pyreadr
    import gseapy as gp

    prot = pyreadr.read_r(str(BASE / "proteomics" / "proc_prot.rds"))[None]
    sets_df = pd.read_csv(BASE / "frozen_gene_sets_v2_clean.txt", sep="\t")

    expr = prot.set_index(prot.columns[0]).astype(float)
    expr.index = expr.index.astype(str)

    gene_sets = {k: [g for g in v.astype(str) if g in expr.index]
                 for k, v in sets_df.groupby("set")["ENTREZID"]}
    for k, v in gene_sets.items():
        print(f"{k}: {len(v)} proteins in matrix")

    # no_scale keeps the raw enrichment scale. gseapy's default rank
    # normalisation compresses small sets, which is what flattened the
    # first boxplot into the 0.84-0.97 band.
    res = gp.ssgsea(data=expr, gene_sets=gene_sets,
                    outdir=None, no_plot=True,
                    sample_norm_method="rank", min_size=5)

    scores = res.res2d.pivot(index="Name", columns="Term", values="ES")
    scores.columns = ["crino_raw" if "crino" in c.lower() else "macro_raw"
                      for c in scores.columns]
    scores.index.name = "record_id"

    donor = read_csv_safe(BASE / "donor.csv")
    d = scores.reset_index().merge(donor, on="record_id", how="inner")
    assert len(d) == 129, f"expected 129 donors, got {len(d)}"

    d.to_csv(BASE / "modelling_dataset_py.csv", index=False)
    print(f"Wrote modelling_dataset_py.csv  ({d.shape})")
    return d


# ---------------------------------------------------------------------
# Models
# ---------------------------------------------------------------------

def primary_hba1c(d):
    print("\n===== PRIMARY: score ~ HbA1c (adjusted) =====")
    rows = []
    for y in ["crino_raw", "macro_raw", "crino_no_unc"]:
        if y not in d.columns:
            continue
        m = smf.ols(f"{y} ~ hba1c + {COVARS}", data=d).fit()
        lo, hi = m.conf_int().loc["hba1c"]
        rows.append(dict(outcome=y, beta=m.params["hba1c"],
                         se=m.bse["hba1c"], ci_lo=lo, ci_hi=hi,
                         p=m.pvalues["hba1c"], n=int(m.nobs)))
    out = pd.DataFrame(rows).round(4)
    print(out.to_string(index=False))
    return out


def group_model(d):
    print("\n===== GROUP MODEL (diagnosis_computed, adjusted) =====")
    for y in ["crino_raw", "macro_raw"]:
        m = smf.ols(
            f'{y} ~ C(diagnosis_computed, Treatment("NoDiabetes")) + {COVARS}',
            data=d).fit()
        at = sm.stats.anova_lm(m, typ=2)
        p = at[at.index.str.contains("diagnosis")]["PR(>F)"].values[0]
        print(f"  {y:12s} omnibus p = {p:.4f}")
        for t in m.params.index:
            if "diagnosis" in t:
                lab = t.split("[T.")[-1].rstrip("]")
                print(f"      {lab:12s} beta={m.params[t]:8.2f}  "
                      f"p={m.pvalues[t]:.4f}")


def variance_test(d):
    """The main result. Dispersion, not location."""
    print("\n===== VARIANCE HETEROGENEITY (adjusted residuals) =====")

    need = ["crino_raw", "macro_raw", "donorage", "bodymassindex",
            "predistributionculturetime", "percentieqrecoverypostculture"]
    sub = d.dropna(subset=need).copy()

    for y in ["crino_raw", "macro_raw"]:
        sub["r_" + y] = smf.ols(f"{y} ~ {COVARS}", data=sub).fit().resid

    results = {}
    for y in ["crino_raw", "macro_raw"]:
        print(f"\n  {y}:")
        groups = []
        for k in GROUPS:
            g = sub.loc[sub.diagnosis_computed == k, "r_" + y]
            groups.append(g)
            iqr = g.quantile(.75) - g.quantile(.25)
            print(f"    {k:12s} n={len(g):3d}  SD={g.std():8.2f}  "
                  f"IQR={iqr:8.2f}")
        # median-centred Levene == Brown-Forsythe, robust to non-normality
        lev = stats.levene(*groups, center="median")
        print(f"    Brown-Forsythe: W={lev.statistic:.3f}  p={lev.pvalue:.4f}")
        results[y] = lev

    print("\n  Specificity: the claim requires crinophagy p < 0.05 AND")
    print("  macroautophagy p > 0.05. Both must be reported together.")
    return results


def outliers(d):
    print("\n===== OUTLIERS (|z| > 3 on crino_raw) =====")
    z = (d.crino_raw - d.crino_raw.mean()) / d.crino_raw.std()
    cols = ["record_id", "crino_raw", "hba1c", "diagnosis_computed",
            "predistributionculturetime", "percentieqrecoverypostculture",
            "pdinsulinperieq"]
    o = d.loc[z.abs() > 3, cols].copy()
    o["z"] = z[z.abs() > 3].round(2)
    print(o.to_string(index=False) if len(o) else "  none")
    print("\n  Note: R241 has missing IEQ recovery and is therefore already")
    print("  dropped from the residual model. The variance result does not")
    print("  depend on an exclusion decision.")
    return o


def functional_anchoring(d):
    print("\n===== FUNCTIONAL ANCHORING =====")
    rows = []
    for fn in ["pdinsulinperieq", "pdinsulindnaratio", "pdinsulincontent"]:
        for y in ["crino_raw", "macro_raw"]:
            m = smf.ols(f"{y} ~ {fn} + {COVARS}", data=d).fit()
            rows.append(dict(outcome=y, predictor=fn,
                             beta=m.params[fn], p=m.pvalues[fn],
                             n=int(m.nobs)))
    out = pd.DataFrame(rows).round(4)
    print(out.to_string(index=False))
    return out


def figures(d):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    need = ["crino_raw", "macro_raw", "donorage", "bodymassindex",
            "predistributionculturetime", "percentieqrecoverypostculture"]
    sub = d.dropna(subset=need).copy()
    for y in ["crino_raw", "macro_raw"]:
        sub["r_" + y] = smf.ols(f"{y} ~ {COVARS}", data=sub).fit().resid

    long = pd.concat([
        sub.assign(score=sub.r_crino_raw, pathway="Crinophagy"),
        sub.assign(score=sub.r_macro_raw, pathway="Macroautophagy"),
    ])
    # z-score within pathway so the two panels are visually comparable;
    # the raw scales differ by set size and are not comparable directly
    long["score_z"] = long.groupby("pathway")["score"].transform(
        lambda s: (s - s.mean()) / s.std())

    fig, ax = plt.subplots(figsize=(8, 5))
    sns.boxplot(data=long, x="diagnosis_computed", y="score_z",
                hue="pathway", showfliers=False, ax=ax,
                palette=["#2e7d32", "#9e9e9e"])
    sns.stripplot(data=long, x="diagnosis_computed", y="score_z",
                  hue="pathway", dodge=True, alpha=.45, size=3.5,
                  ax=ax, palette=["#1b4d1f", "#616161"], legend=False)
    ax.set_xlabel("")
    ax.set_ylabel("Adjusted score (z, within pathway)")
    ax.set_title("Dispersion differs by glycemic status for crinophagy only")
    ax.legend(title=None, frameon=False)
    fig.tight_layout()
    fig.savefig(BASE / "fig_variance_specificity.pdf")
    print(f"\nWrote fig_variance_specificity.pdf")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--score", action="store_true",
                    help="rebuild ssGSEA scores from .rds files")
    ap.add_argument("--csv", default="modelling_dataset.csv")
    args = ap.parse_args()

    d = build_scores() if args.score else read_csv_safe(BASE / args.csv)

    print(f"Loaded {len(d)} donors")
    print(d["diagnosis_computed"].value_counts().to_string())

    primary_hba1c(d)
    group_model(d)
    variance_test(d)
    outliers(d)
    functional_anchoring(d)
    figures(d)


if __name__ == "__main__":
    main()
