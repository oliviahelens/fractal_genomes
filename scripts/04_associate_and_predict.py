"""Genome-wide association of gene-family copy number with leaf dissection, and a
leave-one-clade-out classifier for "fractal" (dissection_order >= 2) leaves.

Association: PGLS (Brownian motion, Grafen branch lengths) per gene family,
    order ~ log2(total genes) + log2(copies + 1)
Prediction: within each training fold, rank families by PGLS on training species
only, take the top K, fit a regularized logistic regression, predict the
held-out clade. Compared with (a) a phylogenetic-kernel baseline and
(b) null traits: label permutations and Brownian-threshold traits simulated on
the same tree.
"""
import pathlib, sys
import numpy as np, pandas as pd, dendropy
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC, RES = ROOT / "data/processed", ROOT / "results"
RES.mkdir(exist_ok=True)
K = int(sys.argv[1]) if len(sys.argv) > 1 else 50
N_NULL = int(sys.argv[2]) if len(sys.argv) > 2 else 100
rng = np.random.default_rng(1)

# ---------- data ----------
ph = pd.read_csv(ROOT / "data/phenotypes/leaf_dissection_curated.tsv", sep="\t", comment="#")
ph = ph[ph.include == 1].set_index("sp")
mat = pd.read_csv(PROC / "family_copy_number.tsv.gz", sep="\t", index_col=0)
total = mat.pop("_total_genes")
species = [s for s in ph.index if s in mat.index]
ph, mat, total = ph.loc[species], mat.loc[species], total.loc[species]
present = (mat > 0).sum()
MIN_PRESENT = 10  # rare families give high-leverage, unstable PGLS fits
mat = mat.loc[:, (present >= MIN_PRESENT) & (mat.std() > 0)]
X = np.log2(mat.values + 1.0)
logtot = np.log2(total.values)
y = ph.dissection_order.values.astype(float)
ybin = (y >= 2).astype(int)
clade = ph.clade.values
n, p = X.shape
print(f"{n} species, {p} families, {ybin.sum()} fractal species in {len(set(clade[ybin == 1]))} clades")

# ---------- phylogenetic covariance (Grafen branch lengths) ----------
tree = dendropy.Tree.get(path=PROC / "species_tree.nwk", schema="newick", preserve_underscores=True)
tree.retain_taxa_with_labels(species)
tree.suppress_unifurcations()
N = len(tree.leaf_nodes())
height = {}
for nd in tree.postorder_node_iter():
    height[nd] = 0.0 if nd.is_leaf() else (len(nd.leaf_nodes()) - 1) / (N - 1)
idx = {s: i for i, s in enumerate(species)}
C = np.eye(n)
for nd in tree.internal_nodes():
    kids = [[idx[l.taxon.label] for l in c.leaf_nodes()] for c in nd.child_nodes()]
    shared = 1.0 - height[nd]
    for a in range(len(kids)):
        for b in range(a + 1, len(kids)):
            C[np.ix_(kids[a], kids[b])] = shared
            C[np.ix_(kids[b], kids[a])] = shared
D = 2.0 - 2.0 * C  # patristic distance on the ultrametric tree


def pgls(Xf, yv, covs, Cm):
    """Vectorized PGLS: t-stat and slope of each column of Xf, controlling for covs."""
    L = np.linalg.cholesky(Cm)
    Wy = np.linalg.solve(L, yv)
    WZ = np.linalg.solve(L, np.column_stack([np.ones(len(yv))] + covs))
    WX = np.linalg.solve(L, Xf)
    Q, _ = np.linalg.qr(WZ)
    ry = Wy - Q @ (Q.T @ Wy)
    rX = WX - Q @ (Q.T @ WX)
    sxx = (rX ** 2).sum(0)
    with np.errstate(divide="ignore", invalid="ignore"):
        b = (rX * ry[:, None]).sum(0) / sxx
        df = len(yv) - WZ.shape[1] - 1
        rv = ((ry ** 2).sum() - b ** 2 * sxx) / df
        t = b / np.sqrt(rv / sxx)
    t = np.nan_to_num(t)
    return t, b, 2 * stats.t.sf(np.abs(t), df)


# ---------- 1. genome-wide association ----------
t, b, pv = pgls(X, y, [logtot], C)
t_ols, _, pv_ols = pgls(X, y, [logtot], np.eye(n))
q = multipletests(pv, method="fdr_bh")[1]
assoc = pd.DataFrame({"family": mat.columns, "t_pgls": t, "slope": b, "p_pgls": pv, "q_pgls": q,
                      "p_ols_no_phylo": pv_ols, "n_species_present": (mat > 0).sum().values,
                      "mean_copies_fractal": mat.values[ybin == 1].mean(0),
                      "mean_copies_other": mat.values[ybin == 0].mean(0)}).sort_values("p_pgls")
assoc["rank"] = np.arange(1, len(assoc) + 1)
assoc.to_csv(RES / "association_pgls.tsv", sep="\t", index=False)
lam = np.median(stats.chi2.isf(pv, 1)) / stats.chi2.ppf(0.5, 1)
lam_ols = np.median(stats.chi2.isf(pv_ols, 1)) / stats.chi2.ppf(0.5, 1)

# positive controls: families of known leaf-dissection regulators (Arabidopsis IDs)
controls = {"STM (KNOX1)": "AT1G62360", "BP/KNAT1 (KNOX1)": "AT4G08150", "KNAT2/6 (KNOX1)": "AT1G70510",
            "CUC2 (NAM)": "AT5G53950", "CUC3": "AT1G76420", "PIN1": "AT1G73590",
            "LMI1 (RCO family, HD-Zip I)": "AT5G03790", "LFY": "AT5G61850", "UFO": "AT1G30950",
            "TFL1": "AT5G03840", "AS1 (ARP)": "AT2G37630", "AS2 (LBD)": "AT1G65620",
            "BOP1": "AT3G57130", "TCP4 (CIN-TCP)": "AT3G15030", "SAW1 (BELL)": "AT1G23380",
            "KNOX2/KNAT3": "AT5G25220", "MIR164-target CUC1": "AT3G15170", "TIE1": "AT4G28840"}
import gzip
g2f = {}
with gzip.open(ROOT / "data/raw/plaza/dicots_HOMFAM.csv.gz", "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        gf, sp, gene = line.rstrip("\n").split("\t")[:3]
        if sp == "ath":
            g2f[gene] = gf
ctrl = []
for name, gene in controls.items():
    fam = g2f.get(gene)
    row = assoc[assoc.family == fam]
    ctrl.append({"gene": name, "ath_id": gene, "family": fam,
                 "rank": int(row["rank"].iloc[0]) if len(row) else None,
                 "percentile": round(100 * row["rank"].iloc[0] / len(assoc), 1) if len(row) else None,
                 "t_pgls": round(float(row.t_pgls.iloc[0]), 2) if len(row) else None,
                 "p_pgls": float(row.p_pgls.iloc[0]) if len(row) else None})
ctrl = pd.DataFrame(ctrl)
ctrl.to_csv(RES / "positive_controls.tsv", sep="\t", index=False)


# ---------- 2. leave-one-clade-out prediction ----------
def cv_predict(yc, Xm=X, k=K):
    yb = (yc >= 2).astype(int)
    score = np.zeros(n)
    for c in np.unique(clade):
        te = clade == c
        tr = ~te
        if yb[tr].sum() == 0:  # nothing to learn from: predict negative
            score[te] = 0.0; continue
        tt, _, _ = pgls(Xm[tr], yc[tr], [logtot[tr]], C[np.ix_(tr, tr)])
        top = np.argsort(-np.abs(tt))[:k]
        F = np.column_stack([Xm[:, top], logtot])
        mu, sd = F[tr].mean(0), F[tr].std(0) + 1e-9
        clf = LogisticRegression(C=0.1, class_weight="balanced", max_iter=2000)
        clf.fit((F[tr] - mu) / sd, yb[tr])
        score[te] = clf.predict_proba((F[te] - mu) / sd)[:, 1]  # probabilities: comparable across folds
    return score


def phylo_kernel(yc):
    yb = (yc >= 2).astype(float)
    score = np.zeros(n)
    for c in np.unique(clade):
        te, tr = clade == c, clade != c
        w = C[np.ix_(te, tr)] + 1e-6
        score[te] = (w @ yb[tr]) / w.sum(1)
    return score


def metrics(yb, s):
    return roc_auc_score(yb, s), average_precision_score(yb, s)


score = cv_predict(y)
auc, ap = metrics(ybin, score)
kauc, kap = metrics(ybin, phylo_kernel(y))
pred = pd.DataFrame({"sp": species, "clade": clade, "dissection_order": y, "fractal": ybin,
                     "cv_score": score}).sort_values("cv_score", ascending=False)
pred["rank"] = np.arange(1, n + 1)
pred.to_csv(RES / "cv_predictions.tsv", sep="\t", index=False)

# null 1: permuted labels; null 2: Brownian-threshold traits with the same number of positives
null_perm, null_bm, assoc_null = [], [], []
L = np.linalg.cholesky(C)
npos3, npos2 = (y == 3).sum(), (y >= 2).sum()
for i in range(N_NULL):
    yp = rng.permutation(y)
    null_perm.append(metrics((yp >= 2).astype(int), cv_predict(yp))[0])
    z = L @ rng.standard_normal(n)
    o = np.argsort(-z)
    ys = np.zeros(n); ys[o[:npos2]] = 2; ys[o[:npos3]] = 3
    ys[o[npos2:npos2 + (y == 1).sum()]] = 1
    null_bm.append(metrics((ys >= 2).astype(int), cv_predict(ys))[0])
    if i < 20:  # calibration of the association test itself
        assoc_null.append(((multipletests(pgls(X, ys, [logtot], C)[2], method="fdr_bh")[1] < 0.05).sum(),
                           (multipletests(pgls(X, yp, [logtot], C)[2], method="fdr_bh")[1] < 0.05).sum()))
null_perm, null_bm = np.array(null_perm), np.array(null_bm)
assoc_null = np.array(assoc_null)

lines = [
    f"species={n} families={p} fractal={ybin.sum()} clades_with_fractal={len(set(clade[ybin == 1]))} K={K} nulls={N_NULL}",
    "",
    "== Association (PGLS, order ~ log2 total genes + log2(copies+1)) ==",
    f"genomic inflation lambda: PGLS {lam:.2f} | ordinary regression (no phylogeny) {lam_ols:.2f}",
    f"families with q<0.05: {(assoc.q_pgls < 0.05).sum()} | q<0.2: {(assoc.q_pgls < 0.2).sum()}",
    f"calibration, families with q<0.05 for random traits (20 each): Brownian-threshold median "
    f"{np.median(assoc_null[:, 0]):.0f} (max {assoc_null[:, 0].max()}), permuted median {np.median(assoc_null[:, 1]):.0f} (max {assoc_null[:, 1].max()})",
    f"families with q<0.05 if phylogeny ignored: {(multipletests(pv_ols, method='fdr_bh')[1] < 0.05).sum()}",
    "top 15:", assoc.head(15)[["family", "t_pgls", "p_pgls", "q_pgls", "n_species_present",
                               "mean_copies_fractal", "mean_copies_other"]].to_string(index=False),
    "",
    "== Positive controls (known leaf-shape regulators) ==", ctrl.to_string(index=False),
    "",
    "== Leave-one-clade-out prediction of fractal (order>=2) ==",
    f"gene-family model: AUROC {auc:.3f}  avg precision {ap:.3f} (chance AP = {ybin.mean():.3f})",
    f"phylogenetic-kernel baseline: AUROC {kauc:.3f}  avg precision {kap:.3f}",
    f"null, permuted labels: AUROC mean {null_perm.mean():.3f}, 95th pct {np.quantile(null_perm, .95):.3f}; "
    f"P(null >= observed) = {(null_perm >= auc).mean():.3f}",
    f"null, Brownian-threshold traits: AUROC mean {null_bm.mean():.3f}, 95th pct {np.quantile(null_bm, .95):.3f}; "
    f"P(null >= observed) = {(null_bm >= auc).mean():.3f}",
    "ranks of fractal species among all held-out predictions:",
    pred[pred.fractal == 1][["sp", "clade", "dissection_order", "cv_score", "rank"]].to_string(index=False),
    "top 15 predicted:", pred.head(15)[["sp", "clade", "dissection_order", "cv_score"]].to_string(index=False),
]
(RES / "summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
