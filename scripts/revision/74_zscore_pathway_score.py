"""Standardized-abundance (mean per-gene z-score) pathway score.
Same residual-dispersion test as the historical ssGSEA/median analyses.
Reviewer 2, major comment 8. Exploratory robustness check.
Usage: python3 -I zscore_score.py <Human Islets dir> <outdir>
"""
import sys, os
import numpy as np, pandas as pd, pyreadr
from scipy.stats import levene

D, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
NBOOT, NPERM, SEED = 5000, 10000, 1
MIN_FRAC = 0.7

def resid(y, X):
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    return y - X @ b

def design(df, cols):
    X = pd.get_dummies(df[cols], drop_first=True, dtype=float)
    return np.column_stack([np.ones(len(df)), X.values])

def zmat(m, robust):
    if robust:
        med = m.median(axis=1); mad = (m.sub(med, axis=0)).abs().median(axis=1) * 1.4826
        return m.sub(med, axis=0).div(mad, axis=0)
    return m.sub(m.mean(axis=1), axis=0).div(m.std(axis=1, ddof=1), axis=0)

def score(m, genes, robust=False):
    g = [x for x in genes if x in m.index]
    z = zmat(m.loc[g], robust)
    ok = z.notna().mean(axis=0) >= MIN_FRAC
    s = z.mean(axis=0, skipna=True); s[~ok] = np.nan
    return s, len(g)

def sdr(r, lab):
    return np.std(r[lab == "T2D"], ddof=1) / np.std(r[lab == "ND"], ddof=1)

def analyse(cohort, setname, variant, s, ngenes, df, cov, status, strata, rng):
    d = df.copy(); d["s"] = s.reindex(d.index).values
    d = d[d["s"].notna()]
    X = design(d, cov); r = resid(d["s"].values, X)
    lab = d[status].values
    two = np.isin(lab, ["ND", "T2D"])
    ratio = sdr(r, lab)
    p_bf2 = levene(r[lab == "ND"], r[lab == "T2D"], center="median").pvalue
    groups = [r[lab == g] for g in pd.unique(lab)]
    p_bf_omni = levene(*groups, center="median").pvalue if len(groups) > 2 else np.nan
    # label permutation (residuals are group-free, so they do not change)
    st = d[strata].values if strata else np.zeros(len(d))
    cnt = 0
    for _ in range(NPERM):
        lp = lab.copy()
        for u in np.unique(st):
            idx = np.where((st == u) & two)[0]
            lp[idx] = rng.permutation(lp[idx])
        if sdr(r, lp) >= ratio: cnt += 1
    p_perm = (cnt + 1) / (NPERM + 1)
    # donor bootstrap, stratified by status (and batch), nuisance refitted
    keys = pd.Series(list(zip(lab, st)))
    bs = []
    for _ in range(NBOOT):
        idx = np.concatenate([rng.choice(np.where(keys == k)[0], (keys == k).sum(), replace=True)
                              for k in keys.unique()])
        db = d.iloc[idx]
        try:
            Xb = design(db, cov)
            if np.linalg.matrix_rank(Xb) < Xb.shape[1]: continue
            rb = resid(db["s"].values, Xb); bs.append(sdr(rb, db[status].values))
        except Exception:
            continue
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return dict(cohort=cohort, set=setname, variant=variant, n_genes=ngenes,
                n_ND=int((lab == "ND").sum()), n_T2D=int((lab == "T2D").sum()), n_fit=len(d),
                sd_ratio=ratio, ci_lo=lo, ci_hi=hi, n_boot_valid=len(bs),
                bf_p_two_group=p_bf2, bf_p_omnibus=p_bf_omni, label_perm_p=p_perm)

rng = np.random.default_rng(SEED)
fs = pd.read_csv(os.path.join(D, "frozen_sets_entrez_to_ensembl.csv"))
gs = {k: v for k, v in fs.groupby("set")}
rows = []

# ---------------- Discovery: processed proteomics -----------------
p = pyreadr.read_r(os.path.join(D, "proteomics/proc_prot.rds"))[None]
p.index = p.iloc[:, 0].astype(str); p = p.iloc[:, 1:].astype(float)
meta = pd.read_csv(os.path.join(D, "donor.csv")).set_index("record_id")
meta["diagnosis_computed"] = meta["diagnosis_computed"].fillna("None")  # pandas reads the literal "None" as NaN
cov = ["donorage", "donorsex", "bodymassindex", "predistributionculturetime", "percentieqrecoverypostculture"]
meta = meta.loc[meta.index.intersection(p.columns)]
meta = meta.dropna(subset=cov + ["diagnosis_computed"])
meta["status"] = meta["diagnosis_computed"].map({"None": "ND", "Pre.T2D": "Pre", "Type2": "T2D"})
p = p[meta.index]
print("discovery donors", len(meta), meta.status.value_counts().to_dict())
hist = pd.read_csv(os.path.join(D, "ssgsea_scores_with_metadata.csv")).set_index("record_id")
for col in ["crinophagy", "macroautophagy"]:
    dd = meta.copy(); dd["s"] = hist[col].reindex(dd.index)
    rr = resid(dd["s"].values, design(dd, cov)); lab = dd.status.values
    print("CHECK ssGSEA", col, "SDR", round(sdr(rr, lab), 3), "BF omnibus",
          round(levene(*[rr[lab == g] for g in pd.unique(lab)], center="median").pvalue, 4))
for setname, key in [("crinophagy", "crinophagy"), ("comparator", "macroautophagy")]:
    ent = gs[key]["ENTREZID"].astype(str).unique().tolist()
    variants = [("z", ent, False), ("robust_z", ent, True)]
    if setname == "crinophagy":
        ins = fs.loc[fs.SYMBOL == "INS", "ENTREZID"].astype(str).tolist()
        variants.append(("z_noINS", [e for e in ent if e not in ins], False))
    for vname, g, rob in variants:
        s, n = score(p, g, rob)
        rows.append(analyse("discovery_processed", setname, vname, s, n, meta, cov, "status", None, rng))
        print(rows[-1])

# ---------------- Replication: GSE164416 VST ----------------------
r = pyreadr.read_r(os.path.join(D, "session_full.RData"))
mat, pd2 = r["mat"], r["pd2"]
pd2 = pd2.copy(); pd2["batch"] = pd2["batch"].astype(str)
mat = mat[pd2.index]
crino = r["crino"].iloc[:, 0].tolist(); macro = r["macro"].iloc[:, 0].tolist()
ins_ens = fs.loc[fs.SYMBOL == "INS", "ENSEMBL"].tolist()
for model, c in [("M1_batch", ["batch"]), ("M2_batch_imm_beta", ["batch", "imm", "beta"])]:
    for setname, ens in [("crinophagy", crino), ("comparator", macro)]:
        variants = [("z", ens, False), ("robust_z", ens, True)]
        if setname == "crinophagy":
            variants.append(("z_noINS", [e for e in ens if e not in ins_ens], False))
        for vname, g, rob in variants:
            s, n = score(mat, g, rob)
            res = analyse("replication_" + model, setname, vname, s, n, pd2, c, "group", "batch", rng)
            rows.append(res); print(res)

out = pd.DataFrame(rows)
out.to_csv(os.path.join(OUT, "74_zscore_pathway_score_results.csv"), index=False)
print(out.round(4).to_string())
