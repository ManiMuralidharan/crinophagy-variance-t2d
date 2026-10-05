#!/usr/bin/env python3
"""
UNC13 paralogue check.

WHAT THIS TESTS. The four Munc13 paralogues are not interchangeable:
  UNC13A/B/C (Munc13-1/-2/-3) are C1-domain, DAG-regulated priming
    factors at neuronal and neuroendocrine active zones
  UNC13D    (Munc13-4) lacks the C1 domain and acts in lysosome-related
    organelle fusion, which is the crinophagy-relevant one

So a variance signal in UNC13A would be a DIFFERENT finding, not a
rescue of the UNC13D result. The real value here is as a technical
control: if all paralogues show the same pattern, it points to
family-level detectability rather than Munc13-4-specific biology.

Usage:
  cd "/Users/mmani3/Downloads/Human Islets"
  python3 07_unc13_paralogues.py
"""

from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats

BASE = Path("/Users/mmani3/Downloads/Human Islets")

COVARS = ("donorage + C(donorsex) + bodymassindex "
          "+ predistributionculturetime + percentieqrecoverypostculture")

GROUPS = ["NoDiabetes", "Pre.T2D", "Type2"]

NEED = ["crino_raw", "donorage", "bodymassindex",
        "predistributionculturetime", "percentieqrecoverypostculture"]

NUMERIC = ["crino_raw", "macro_raw", "hba1c", "donorage", "bodymassindex",
           "predistributionculturetime", "percentieqrecoverypostculture"]

# Entrez and canonical UniProt for each paralogue, plus RAB3A as the
# positive reference (clean: 129/129 detection, SD ratio 1.73)
PARALOGUES = [
    ("UNC13A", "23025",  "Q9UPW8"),
    ("UNC13B", "10497",  "O14795"),
    ("UNC13C", "440279", "Q8NB66"),
    ("UNC13D", "201294", "Q70J99"),
    ("RAB3A",  "5864",   "P20336"),
]


def read_csv_safe(path):
    """pandas reads 'None' as NaN, dropping all 91 non-diabetic donors."""
    d = pd.read_csv(path, na_values=["NA", ""], keep_default_na=False)
    for c in NUMERIC:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d["diagnosis_computed"] = d["diagnosis_computed"].replace(
        "None", "NoDiabetes")
    return d


def main():
    import pyreadr

    d = read_csv_safe(BASE / "modelling_dataset.csv")

    proc = pyreadr.read_r(str(BASE / "proteomics" / "proc_prot.rds"))[None]
    proc = proc.set_index(proc.columns[0])
    proc.index = proc.index.astype(str)
    proc = proc.astype(float)

    raw = pyreadr.read_r(str(BASE / "proteomics" / "unproc_prot.rds"))[None]
    raw = raw.set_index(raw.columns[0])
    raw.index = raw.index.astype(str)
    raw = raw.astype(float)

    # ---- 1. Detection ------------------------------------------------
    print("=" * 64)
    print("1. DETECTION IN THE RAW (UniProt-keyed) MATRIX")
    print("=" * 64)

    grp = d.set_index("record_id")["diagnosis_computed"]

    for sym, eid, uni in PARALOGUES:
        if uni not in raw.index:
            print(f"\n  {sym:8s} accession {uni} ABSENT from raw matrix")
            print(f"           (not quantified in this experiment)")
            continue

        v = raw.loc[uni]
        det = (~v.isna()) & (v != 0)
        print(f"\n  {sym:8s} {uni}   detected in "
              f"{det.sum()}/{len(v)} ({100*det.mean():.1f}%)")

        # detection by group: this is what compromised UNC13D
        ids = [str(i) for i in v.index]
        m = pd.DataFrame({"record_id": ids, "detected": det.to_numpy()})
        m = m[m.record_id.isin(grp.index)]
        m["group"] = grp.loc[m.record_id].to_numpy()

        tab = pd.crosstab(m.group, m.detected).reindex(GROUPS).fillna(0)
        for g in GROUPS:
            r = m[m.group == g]
            if len(r):
                print(f"           {g:12s} {r.detected.sum():3.0f}/{len(r):3d}"
                      f" = {100*r.detected.mean():5.1f}%")

        if tab.shape[1] > 1 and (tab.to_numpy() > 0).all():
            chi = stats.chi2_contingency(tab.to_numpy())
            flag = "  <-- CONFOUNDED" if chi[1] < 0.05 else ""
            print(f"           chi-square p = {chi[1]:.4f}{flag}")
        else:
            print(f"           (no group variation in detection)")

    # ---- 2. Variance -------------------------------------------------
    print("\n" + "=" * 64)
    print("2. VARIANCE HETEROGENEITY (processed matrix, adjusted)")
    print("=" * 64)

    sub = d.dropna(subset=NEED).copy()
    keep = [r for r in sub.record_id if r in proc.columns]
    meta = sub.set_index("record_id").loc[keep]
    rng = np.random.default_rng(1)

    rows = []
    for sym, eid, uni in PARALOGUES:
        if eid not in proc.index:
            print(f"\n  {sym:8s} not in processed matrix "
                  f"(filtered out or never measured)")
            continue

        df = meta.assign(y=proc.loc[eid, keep].to_numpy()).dropna(subset=["y"])
        if df.y.std() == 0 or len(df) < 20:
            print(f"\n  {sym:8s} insufficient data (n={len(df)})")
            continue

        df["r"] = smf.ols(f"y ~ {COVARS}", data=df).fit().resid
        a = df.loc[df.diagnosis_computed == "Type2", "r"]
        b = df.loc[df.diagnosis_computed == "NoDiabetes", "r"]
        if len(a) < 3 or b.std() == 0:
            continue

        ratio = a.std(ddof=1) / b.std(ddof=1)
        p = stats.levene(a, b, center="median").pvalue

        boot = np.array([
            rng.choice(a, len(a), True).std(ddof=1) /
            rng.choice(b, len(b), True).std(ddof=1) for _ in range(5000)])
        lo, hi = np.percentile(boot, [2.5, 97.5])

        rows.append(dict(symbol=sym, n=len(df), sd_ratio=ratio,
                         ci_lo=lo, ci_hi=hi, p=p))

    out = pd.DataFrame(rows)
    print()
    print(out.round(4).to_string(index=False))

    # ---- 3. Interpretation -------------------------------------------
    print("\n" + "=" * 64)
    print("3. HOW TO READ THIS")
    print("=" * 64)

    fam = out[out.symbol.str.startswith("UNC13")]
    if len(fam) > 1:
        n_elev = (fam.sd_ratio > 1.3).sum()
        print(f"\n  Paralogues with SD ratio > 1.3: {n_elev}/{len(fam)}")
        if n_elev >= 3:
            print("  => family-wide pattern. Points to shared detectability")
            print("     or expression behaviour, NOT Munc13-4-specific")
            print("     biology. Weakens the UNC13D claim further.")
        elif n_elev == 1:
            print("  => isolated to one paralogue. More consistent with a")
            print("     specific effect, though the UNC13D detection")
            print("     confound (chi-square p = 0.015) still stands.")

    print("\n  Remember: UNC13A/B/C are Munc13-1/-2/-3, priming factors at")
    print("  the active zone. They are NOT substitutes for Munc13-4. A hit")
    print("  in UNC13A is a separate finding about granule priming, and")
    print("  would need its own framing and its own comparator.")
    print("\n  RAB3A remains the clean result: 129/129 detection, no")
    print("  dropout confound, bootstrap CI excluding 1.0.")

    out.to_csv(BASE / "unc13_paralogue_variance.csv", index=False)
    print("\n  Wrote unc13_paralogue_variance.csv")


if __name__ == "__main__":
    main()
