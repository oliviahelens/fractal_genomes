"""Relative evolutionary rates (RER, after RERconverge) vs leaf dissection.

For each gene tree (branch lengths on the fixed species topology), each branch
length is compared with the genome-wide expectation for that branch (master
tree, from a concatenated alignment): RER = residual of sqrt(gene length) on
sqrt(master length), scaled for heteroscedasticity. Terminal-branch RERs are
then centered per species (removes genome-wide rate effects such as life
history or assembly quality).

Tests, per gene: are fractal (order >= 2) lineages faster/slower (AUC), and does
RER track dissection order (Spearman)? Calibrated by "permulations": Brownian-
motion traits simulated on the species tree, thresholded to the same number of
fractal tips. Then a leave-one-clade-out classifier on terminal RERs.

Usage: python3 07_rer.py [n_null_assoc] [n_null_cv] [K]
"""
import pathlib, sys, subprocess, gzip
import numpy as np, pandas as pd, dendropy
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC, RES = ROOT / "data/processed", ROOT / "results/rer"
WORK = ROOT / "data/work/gene_trees"
RES.mkdir(parents=True, exist_ok=True)
BIN = pathlib.Path("/home/user/tools/env/bin")
N_NULL = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
N_NULL_CV = int(sys.argv[2]) if len(sys.argv) > 2 else 100
K = int(sys.argv[3]) if len(sys.argv) > 3 else 50
rng = np.random.default_rng(7)

ph = pd.read_csv(ROOT / "data/phenotypes/leaf_dissection_curated.tsv", sep="\t", comment="#")
ph = ph[ph.include == 1].set_index("sp")
species = ph.index.tolist(); n = len(species); idx = {s: i for i, s in enumerate(species)}
y = ph.dissection_order.values.astype(float); ybin = (y >= 2).astype(int); clade = ph.clade.values

tree = dendropy.Tree.get(path=PROC / "species_tree.nwk", schema="newick", preserve_underscores=True)
tree.retain_taxa_with_labels(species); tree.suppress_unifurcations()
ALL = frozenset(species)


def canon(tips, universe):
    """Unrooted split key: the side not containing the universe's alphabetically first tip."""
    ref = min(universe)
    return frozenset(tips) if ref not in tips else frozenset(universe - set(tips))


full_splits = []  # tip sets of every edge in the rooted species tree
for nd in tree.preorder_node_iter():
    if nd.parent_node is not None:
        full_splits.append(frozenset(l.taxon.label for l in nd.leaf_nodes()))


def tree_lengths(path):
    t = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True)
    tips = frozenset(l.taxon.label for l in t.leaf_node_iter())
    out = {}
    for nd in t.preorder_node_iter():
        if nd.parent_node is None:
            continue
        k = canon({l.taxon.label for l in nd.leaf_nodes()}, tips)
        if 0 < len(k) < len(tips):
            out[k] = out.get(k, 0.0) + (nd.edge.length or 0.0)
    return tips, out


done = sorted(p.parent.name for p in WORK.glob("*/done"))
print(f"{len(done)} gene trees", flush=True)

# ---------- master tree: concatenation of the 300 most complete genes ----------
master_path = RES / "master.treefile"
cov = {f: sum(1 for l in open(WORK / f / "rep.trim.aln") if l.startswith(">")) for f in done}
if not master_path.exists():
    top = sorted(done, key=lambda f: -cov[f])[:300]
    seqs = {s: [] for s in species}
    for f in top:
        aln, name = {}, None
        for line in open(WORK / f / "rep.trim.aln"):
            if line.startswith(">"): name = line[1:].strip(); aln[name] = []
            else: aln[name].append(line.strip())
        L = len("".join(next(iter(aln.values()))))
        for s in species:
            seqs[s].append("".join(aln[s]) if s in aln else "-" * L)
    cat = RES / "master_concat.fa"
    with open(cat, "w") as fh:
        for s in species: fh.write(f">{s}\n{''.join(seqs[s])}\n")
    (RES / "master_topology.nwk").write_text(tree.as_string(schema="newick").replace("[&R] ", ""))
    subprocess.run([BIN / "iqtree", "-s", cat, "-te", RES / "master_topology.nwk", "-m", "LG+G4", "-nt", "4",
                    "-quiet", "-redo", "-pre", RES / "master"], check=True, capture_output=True)
    cat.unlink()
_, master_full = tree_lengths(master_path)
master_by_edge = {e: master_full[canon(e, ALL)] for e in full_splits}  # rooted edges -> length (root edge shared)
# the two root edges map to one unrooted split; give each half so sums stay right
root_kids = [frozenset(l.taxon.label for l in c.leaf_nodes()) for c in tree.seed_node.child_nodes()]
for e in root_kids:
    master_by_edge[e] = master_full[canon(e, ALL)] / len(root_kids)

# ---------- RER per gene ----------
rows_x, rows_y, rows_gene, rows_split = [], [], [], []
for f in done:
    tips, gl = tree_lengths(WORK / f / "iq.treefile")
    mp = {}
    for e, L in master_by_edge.items():
        r = e & tips
        if 0 < len(r) < len(tips):
            k = canon(r, tips); mp[k] = mp.get(k, 0.0) + L
    keys = [k for k in gl if k in mp]
    rows_x += [mp[k] for k in keys]; rows_y += [gl[k] for k in keys]
    rows_gene += [f] * len(keys); rows_split += keys
df = pd.DataFrame({"gene": rows_gene, "x": np.sqrt(rows_x), "y": np.sqrt(rows_y)})
df["split"] = rows_split
res = np.empty(len(df))
for g, sub in df.groupby("gene").indices.items():
    xx, yy = df.x.values[sub], df.y.values[sub]
    A = np.column_stack([np.ones(len(xx)), xx])
    beta = np.linalg.lstsq(A, yy, rcond=None)[0]
    res[sub] = yy - A @ beta
df["res"] = res
df["bin"] = pd.qcut(df.x, 20, labels=False, duplicates="drop")
df["rer"] = df.res / df.groupby("bin").res.transform("std")
term = df[df.split.map(len) == 1].copy()
term["sp"] = term.split.map(lambda s: next(iter(s)))
R = term.pivot_table(index="gene", columns="sp", values="rer").reindex(columns=species)
R = R - R.median(axis=0)  # remove per-species genome-wide rate offset
R.to_csv(RES / "terminal_rer.tsv.gz", sep="\t", compression="gzip")
Rv = R.values; G = Rv.shape[0]
print(f"RER matrix {R.shape}, {np.isnan(Rv).mean():.1%} missing", flush=True)

# ---------- permulation-calibrated gene tests ----------
height = {}
Ntips = len(tree.leaf_nodes())
for nd in tree.postorder_node_iter():
    height[nd] = 0.0 if nd.is_leaf() else (len(nd.leaf_nodes()) - 1) / (Ntips - 1)
C = np.eye(n)
for nd in tree.internal_nodes():
    kids = [[idx[l.taxon.label] for l in c.leaf_nodes()] for c in nd.child_nodes()]
    for a in range(len(kids)):
        for b in range(a + 1, len(kids)):
            C[np.ix_(kids[a], kids[b])] = 1 - height[nd]; C[np.ix_(kids[b], kids[a])] = 1 - height[nd]
Lc = np.linalg.cholesky(C)
mask = ~np.isnan(Rv)
ranks = np.full_like(Rv, np.nan)
for i in range(G):
    ranks[i, mask[i]] = stats.rankdata(Rv[i, mask[i]])
npres = mask.sum(1)


def auc_all(lab):
    """AUC per gene of RER separating lab==1 tips (Mann-Whitney), NaN-aware."""
    pos = mask & (lab[None, :] == 1)
    n1 = pos.sum(1); n0 = npres - n1
    rs = np.where(pos, ranks, 0).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (rs - n1 * (n1 + 1) / 2) / (n1 * n0)


def spear_all(yv):
    out = np.empty(G)
    for i in range(G):
        m = mask[i]; out[i] = stats.spearmanr(Rv[i, m], yv[m])[0]
    return out


def sim_trait():
    z = Lc @ rng.standard_normal(n); o = np.argsort(-z)
    ys = np.zeros(n); ys[o[:(y >= 1).sum()]] = 1; ys[o[:(y >= 2).sum()]] = 2; ys[o[:(y >= 3).sum()]] = 3
    return ys


obs_auc = auc_all(ybin)
obs_rho = spear_all(y)
null_auc = np.empty((N_NULL, G)); null_rho = np.empty((N_NULL, G)); sims = []
for k in range(N_NULL):
    ys = sim_trait(); sims.append(ys)
    null_auc[k] = auc_all((ys >= 2).astype(int))
    if k < 200: null_rho[k] = spear_all(ys)
null_rho = null_rho[:min(N_NULL, 200)]


def calib(obs, null):
    mu, sd = np.nanmean(null, 0), np.nanstd(null, 0)
    z = (obs - mu) / sd
    p_emp = (1 + (np.abs(null - mu) >= np.abs(obs - mu)).sum(0)) / (1 + len(null))
    return z, 2 * stats.norm.sf(np.abs(z)), p_emp


z_auc, p_auc, pe_auc = calib(obs_auc, null_auc)
z_rho, p_rho, pe_rho = calib(obs_rho, null_rho)
out = pd.DataFrame({"family": R.index, "n_species": npres, "auc_fractal": obs_auc, "z_auc": z_auc,
                    "p_auc": p_auc, "p_auc_empirical": pe_auc, "q_auc": multipletests(np.nan_to_num(p_auc, nan=1), method="fdr_bh")[1],
                    "spearman_order": obs_rho, "z_rho": z_rho, "p_rho": p_rho,
                    "q_rho": multipletests(np.nan_to_num(p_rho, nan=1), method="fdr_bh")[1]}).sort_values("p_auc")
# Arabidopsis members for interpretation
ga = pd.read_csv(PROC / "gene_family_assignments.tsv.gz", sep="\t")
ath = ga[(ga.sp == "ath") & ga.gf.isin(R.index)].groupby("gf").gene.apply(lambda s: ",".join(sorted(s)[:4]))
out["ath_genes"] = out.family.map(ath).fillna("")
out.to_csv(RES / "rer_association.tsv", sep="\t", index=False)

# how many hits would random traits give? (same procedure applied to simulated traits)
null_hits = []
for k in range(min(20, N_NULL)):
    zz = (null_auc[k] - np.nanmean(np.delete(null_auc, k, 0), 0)) / np.nanstd(np.delete(null_auc, k, 0), 0)
    null_hits.append((multipletests(np.nan_to_num(2 * stats.norm.sf(np.abs(zz)), nan=1), method="fdr_bh")[1] < 0.1).sum())


# ---------- leave-one-clade-out classifier on terminal RER ----------
Xall = np.nan_to_num(Rv.T)  # species x genes, missing = 0 (= no deviation)


def cv(yc):
    yb = (yc >= 2).astype(int); score = np.zeros(n)
    for c in np.unique(clade):
        te, tr = clade == c, clade != c
        if yb[tr].sum() == 0:
            continue
        lab = np.where(tr, yb, -1)
        pos = mask & tr[None, :] & (yb[None, :] == 1); neg = mask & tr[None, :] & (yb[None, :] == 0)
        # training-only AUC per gene (ranks within training species)
        a = np.empty(G)
        for i in range(G):
            m = mask[i] & tr
            r = stats.rankdata(Rv[i, m]); l = yb[m]
            n1, n0 = l.sum(), (1 - l).sum()
            a[i] = (r[l == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0) if n1 and n0 else 0.5
        top = np.argsort(-np.abs(a - 0.5))[:K]
        F = Xall[:, top]
        clf = LogisticRegression(C=0.1, class_weight="balanced", max_iter=2000).fit(F[tr], yb[tr])
        score[te] = clf.predict_proba(F[te])[:, 1]
    return score


score = cv(y)
auc, ap = roc_auc_score(ybin, score), average_precision_score(ybin, score)
null_cv = np.array([roc_auc_score((s >= 2).astype(int), cv(s)) for s in sims[:N_NULL_CV]])
pred = pd.DataFrame({"sp": species, "clade": clade, "order": y, "fractal": ybin, "score": score}).sort_values("score", ascending=False)
pred["rank"] = np.arange(1, n + 1)
pred.to_csv(RES / "rer_cv_predictions.tsv", sep="\t", index=False)

lines = [
    f"gene trees: {G}; species {n}; fractal {ybin.sum()} in {len(set(clade[ybin == 1]))} clades; nulls {N_NULL} (assoc), {N_NULL_CV} (cv)",
    "",
    "== Genes whose rate differs in fractal lineages (terminal-branch RER, permulation-calibrated) ==",
    f"q<0.1: {(out.q_auc < 0.1).sum()}  (random traits through the same pipeline: median {np.median(null_hits):.0f}, max {max(null_hits)})",
    f"q<0.1 for correlation with dissection order: {(out.q_rho < 0.1).sum()}",
    out.head(20)[["family", "n_species", "auc_fractal", "z_auc", "p_auc", "q_auc", "spearman_order", "ath_genes"]].to_string(index=False),
    "",
    "== Leave-one-clade-out prediction from terminal RER ==",
    f"AUROC {auc:.3f}  avg precision {ap:.3f} (chance {ybin.mean():.3f})",
    f"null (Brownian-threshold traits): AUROC mean {null_cv.mean():.3f}, 95th pct {np.quantile(null_cv, .95):.3f}; P(null >= obs) = {(null_cv >= auc).mean():.3f}",
    pred[pred.fractal == 1][["sp", "clade", "order", "score", "rank"]].to_string(index=False),
]
(RES / "summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
