"""Regenerate Figure 2C (1,000,000-set competitive null, Model 2) and
Supplementary Figure S4 (processed group means, no positivity filter or log10).
Usage: python3 -I regen.py <Human Islets dir> <outdir>
"""
import sys, os
import numpy as np, pandas as pd, pyreadr
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm

D, OUT = sys.argv[1], sys.argv[2]
os.makedirs(OUT, exist_ok=True)
for f in fm.findSystemFonts():
    if "LiberationSans" in f: fm.fontManager.addfont(f)
plt.rcParams.update({"font.family": "Liberation Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
GREEN, GREY, ORANGE = "#1b5e20", "#90a4ae", "#ef6c00"

# ---------------- Figure 2C ----------------
r = pyreadr.read_r(os.path.join(D, "session_gse164416.RData"))
bg = r["bg2"].iloc[:, 0].values
bg = bg[np.isfinite(bg)]
full = pd.read_csv(os.path.join(D, "replication_GSE164416/canonical_full.csv"))
m2 = full[full.model == "M2_full"]
NSIM, rng = 1_000_000, np.random.default_rng(1)

def null_medians(k):
    out = np.empty(NSIM); done = 0
    while done < NSIM:
        n = min(200_000, NSIM - done)
        idx = rng.integers(0, len(bg), size=(n, k))
        s = np.sort(idx, axis=1)
        dup = (np.diff(s, axis=1) == 0).any(axis=1)
        while dup.any():  # redraw rows with repeated genes: sampling without replacement
            idx[dup] = rng.integers(0, len(bg), size=(dup.sum(), k))
            s = np.sort(idx[dup], axis=1)
            sub = (np.diff(s, axis=1) == 0).any(axis=1)
            tmp = np.zeros_like(dup); tmp[np.where(dup)[0][sub]] = True; dup = tmp
        out[done:done + n] = np.median(bg[idx], axis=1); done += n
    return out

rows, nulls = [], {}
for s, lab in [("crinophagy", "Crinophagy set (31 genes, frozen)"),
               ("comparator", "Macroautophagy comparator (14 genes, frozen)")]:
    g = m2.loc[m2.set == s, "sd_ratio"].dropna().values
    obs = np.median(g); nl = null_medians(len(g))
    p = (np.sum(nl >= obs) + 1) / (NSIM + 1)
    rows.append(dict(set=s, n_genes=len(g), observed_median=obs, null_median=np.median(nl),
                     competitive_p=p, n_random_sets=NSIM, seed=1)); nulls[lab] = (nl, obs, p)
res = pd.DataFrame(rows); print(res.to_string())
res.to_csv(os.path.join(OUT, "76_fig2C_competitive_1e6.csv"), index=False)

fig, axs = plt.subplots(1, 2, figsize=(8, 3.3))
for ax, (lab, (nl, obs, p)) in zip(axs, nulls.items()):
    ax.hist(nl, bins=80, color="#dcdcdc", edgecolor="none")
    ax.axvline(obs, color=GREEN, lw=1.4)
    ax.text(0.03, 0.95, f"observed median {obs:.3f}\ncompetitive P = {p:.3f}",
            transform=ax.transAxes, ha="left", va="top", color=GREEN, fontsize=8)
    ax.set_title(lab, fontsize=8.5, fontweight="bold")
    ax.set_xlabel("Median SD ratio of a size-matched random gene set")
    ax.set_ylabel("Random gene sets")
fig.suptitle("C  Competitive random-gene-set comparison (Model 2)", x=0.01, ha="left",
             fontweight="bold", fontsize=10)
fig.text(0.01, -0.02, "1,000,000 size-matched random gene sets drawn without replacement from the Model 2 "
         "genome-wide SD-ratio distribution (20,255 genes); matched for size only.",
         fontsize=7, color="#555555")
fig.tight_layout(rect=(0, 0.02, 1, 0.95))
for ext in ("pdf", "png", "svg"):
    fig.savefig(os.path.join(OUT, f"Figure2C_competitive_1e6.{ext}"), dpi=300, bbox_inches="tight")
plt.close(fig)

# ---------------- Supplementary Figure S4 ----------------
prot = pyreadr.read_r(os.path.join(D, "proteomics/proc_prot.rds"))[None]
prot.index = prot.iloc[:, 0].astype(str); prot = prot.iloc[:, 1:].astype(float)
rna = pyreadr.read_r(os.path.join(D, "Bulk RNA/proc_rnaseq.rds"))[None]
rna.index = rna["gene_id"].astype(str)
rna = rna.drop(columns=[c for c in ["accession", "gene_id", "symbol", "name", "biotype"] if c in rna.columns]).astype(float)
rna = rna[~rna.index.duplicated()]
donor = pd.read_csv(os.path.join(D, "donor.csv")).set_index("record_id")
donor["diagnosis_computed"] = donor["diagnosis_computed"].fillna("None")
sets = pd.read_csv(os.path.join(D, "frozen_sets_entrez_to_ensembl.csv"))
key = sets[["SYMBOL", "ENTREZID", "set"]].drop_duplicates("SYMBOL").copy()
key["set"] = key["set"].replace({"macroautophagy": "comparator"})
key = pd.concat([key, pd.DataFrame([dict(SYMBOL="MAP1LC3B", ENTREZID=81631, set="comparator"), dict(SYMBOL="STX1A", ENTREZID=6804, set="post-hoc")])])
key["ENTREZID"] = key["ENTREZID"].astype(int).astype(str)
print("panel genes:", len(key), key.set.value_counts().to_dict())

shared = [d for d in prot.columns if d in rna.columns and d in donor.index]
grp = donor.loc[shared, "diagnosis_computed"].map({"None": "Non-diabetic", "Pre.T2D": "Pre-T2D", "Type2": "T2D"})
print("donors with both layers:", len(shared), grp.value_counts().to_dict())
both = key[key.ENTREZID.isin(prot.index) & key.ENTREZID.isin(rna.index)]
print("panel genes in both processed layers:", len(both), "of", len(key))
missing = key[~key.SYMBOL.isin(both.SYMBOL)].SYMBOL.tolist()

recs = []
for g in ["Non-diabetic", "Pre-T2D", "T2D"]:
    dd = grp[grp == g].index
    for _, k in both.iterrows():
        recs.append(dict(group=g, n_donors=len(dd), SYMBOL=k.SYMBOL, set=k.set,
                         protein_mean=prot.loc[k.ENTREZID, dd].mean(),
                         rna_mean=rna.loc[k.ENTREZID, dd].mean()))
df = pd.DataFrame(recs).dropna()
df.to_csv(os.path.join(OUT, "77_figS4_group_means.csv"), index=False)
pd.Series(missing, name="not_in_both_processed_layers").to_csv(
    os.path.join(OUT, "77_figS4_genes_not_plotted.csv"), index=False)

col = {"crinophagy": GREEN, "comparator": GREY, "post-hoc": ORANGE}
LABEL = {"INS", "CHGA", "UNC13D", "RAB3A", "RAB27A", "STX1A", "CTSD", "GABARAP"}
fig, axs = plt.subplots(1, 3, figsize=(10, 3.6), sharex=True, sharey=True)
for ax, g in zip(axs, ["Non-diabetic", "Pre-T2D", "T2D"]):
    d = df[df.group == g]
    for s, c in col.items():
        e = d[d.set == s]
        ax.scatter(e.protein_mean, e.rna_mean, s=16, color=c, label=s, zorder=3)
    for _, e in d[d.SYMBOL.isin(LABEL)].iterrows():
        ax.annotate(e.SYMBOL, (e.protein_mean, e.rna_mean), fontsize=6.5,
                    xytext=(3, 3), textcoords="offset points", color="#333333")
    rho = d[["protein_mean", "rna_mean"]].corr(method="spearman").iloc[0, 1]
    ax.set_title(f"{g} (n = {d.n_donors.iloc[0]})", fontsize=9, fontweight="bold")
    ax.text(0.03, 0.97, f"Spearman rho = {rho:.2f}\n{len(d)} genes", transform=ax.transAxes,
            va="top", fontsize=7, color="#555555")
    ax.set_xlabel("Mean processed protein abundance")
axs[0].set_ylabel("Mean processed RNA abundance")
axs[2].legend(frameon=False, fontsize=7, loc="lower right")
fig.suptitle("S4  Transcript and protein abundance across the frozen-set panel", x=0.01,
             ha="left", fontweight="bold", fontsize=10)
fig.tight_layout(rect=(0, 0, 1, 0.94))
for ext in ("pdf", "png", "svg"):
    fig.savefig(os.path.join(OUT, f"FigureS4_group_means.{ext}"), dpi=300, bbox_inches="tight")
print("S4 genes plotted per group:", df.groupby("group").size().to_dict(), "missing:", missing)
