"""Relative evolutionary rates vs leaf dissection on the expanded species set
(PLAZA + NCBI genomes + 1KP transcriptomes). Same method as 07_rer.py, plus:
  - species filter: phenotype include=1, in the expanded tree, >= MIN_COMP of
    the families recovered;
  - clades for leave-one-out = NCBI orders;
  - checks for data-type artifacts: (a) can terminal RER predict genome vs
    transcriptome? (b) per-gene tests repeated on genomes only; (c) repeated
    without low-confidence phenotype scores.
Usage: python3 14_rer_expanded.py [n_null_assoc] [n_null_cv] [K]
"""
import pathlib, sys, subprocess
import numpy as np, pandas as pd, dendropy
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC, RES = ROOT / "data/processed", ROOT / "results/rer_expanded"
WORK = ROOT / "data/work/gene_trees_x"
RES.mkdir(parents=True, exist_ok=True)
BIN = pathlib.Path("/home/user/tools/env/bin")
N_NULL = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
N_NULL_CV = int(sys.argv[2]) if len(sys.argv) > 2 else 100
K = int(sys.argv[3]) if len(sys.argv) > 3 else 50
MIN_COMP, MIN_SP_GENE, MAX_COLS = 0.5, 100, 30000
rng = np.random.default_rng(7)

# ---------- species and phenotypes ----------
st = pd.read_csv(PROC / "species_table.tsv", sep="\t").set_index("sp")
ph = pd.concat([pd.read_csv(ROOT / "data/phenotypes/leaf_dissection_curated.tsv", sep="\t", comment="#")] +
               [pd.read_csv(p, sep="\t") for p in sorted((ROOT / "data/phenotypes/curation_parts").glob("scored_*.tsv"))])
ph = ph.drop_duplicates("sp").set_index("sp")
comp = pd.read_csv(PROC / "completeness.tsv", sep="\t").set_index("sp").completeness
tree = dendropy.Tree.get(path=PROC / "species_tree_expanded_binary.nwk", schema="newick", preserve_underscores=True,
                         rooting="force-rooted")
tips = {l.taxon.label for l in tree.leaf_node_iter()}
done = sorted(p.parent.name for p in WORK.glob("*/done"))
ok = ph[(ph.include == 1) & ph.index.isin(tips)].index
ok = [s for s in ok if st.loc[s, "source"] == "plaza" or comp.get(s, 0) >= MIN_COMP]
species = sorted(ok); n = len(species); idx = {s: i for i, s in enumerate(species)}
y = ph.loc[species, "dissection_order"].values.astype(float); ybin = (y >= 2).astype(int)
order = st.loc[species, "order"].fillna("NA").values
dtype = st.loc[species, "data_type"].values
lowconf = (ph.loc[species, "confidence"] == "low").values
tree.retain_taxa_with_labels(species); tree.suppress_unifurcations()
ALL = frozenset(species)
print(f"{len(done)} gene trees; {n} species ({(dtype == 'genome').sum()} genomes); fractal {ybin.sum()} in "
      f"{len(set(order[ybin == 1]))} orders", flush=True)


def canon(tips_, universe):
    ref = min(universe)
    return frozenset(tips_) if ref not in tips_ else frozenset(universe - set(tips_))


full_splits = [frozenset(l.taxon.label for l in nd.leaf_nodes()) for nd in tree.preorder_node_iter() if nd.parent_node is not None]


def tree_lengths(path, keep=None):
    t = dendropy.Tree.get(path=path, schema="newick", preserve_underscores=True)
    if keep is not None:
        lab = [l.taxon.label for l in t.leaf_node_iter() if l.taxon.label in keep]
        t.retain_taxa_with_labels(lab); t.suppress_unifurcations()
    tp = frozenset(l.taxon.label for l in t.leaf_node_iter())
    out = {}
    for nd in t.preorder_node_iter():
        if nd.parent_node is None: continue
        k = canon({l.taxon.label for l in nd.leaf_nodes()}, tp)
        if 0 < len(k) < len(tp): out[k] = out.get(k, 0.0) + (nd.edge.length or 0.0)
    return tp, out


# ---------- master tree: concatenation of the 300 most complete genes (column-subsampled) ----------
def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"): name = line[1:]; out[name] = []
        else: out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


master_path = RES / "master.nwk"
if not master_path.exists():
    cov = {f: sum(1 for l in open(WORK / f / "aln.fa") if l.startswith(">") and l[1:].strip() in ALL) for f in done}
    top = sorted(done, key=lambda f: -cov[f])[:300]
    seqs = {s: [] for s in species}
    for f in top:
        aln = read_fasta(WORK / f / "aln.fa"); L = len(next(iter(aln.values())))
        for s in species: seqs[s].append(aln.get(s, "-" * L))
    # FastTree memory scales with taxa x columns: keep columns >= 50% filled, then a seeded
    # random sample of MAX_COLS of them (~5 GB for ~1,100 taxa), enough for branch lengths
    M = np.array([list("".join(seqs[s])) for s in species])
    good = np.where((M != "-").mean(0) >= 0.5)[0]
    cols = np.sort(np.random.default_rng(0).choice(good, min(MAX_COLS, len(good)), replace=False))
    cat = RES / "master_concat.fa"
    with open(cat, "w") as fh:
        for i, s in enumerate(species): fh.write(f">{s}\n{''.join(M[i, cols])}\n")
    print(f"master alignment: {len(cols)} of {M.shape[1]} columns ({len(good)} >= 50% filled)", flush=True)
    del M
    (RES / "master_topology.nwk").write_text(tree.as_string(schema="newick").replace("[&R] ", "").strip())
    tmp = master_path.with_suffix(".tmp")
    with open(tmp, "w") as out:
        subprocess.run([BIN / "FastTreeMP", "-lg", "-gamma", "-nome", "-mllen", "-intree", RES / "master_topology.nwk",
                        "-quiet", "-nopr", cat], stdout=out, stderr=subprocess.DEVNULL, check=True)
    tmp.rename(master_path); cat.unlink()
_, master_full = tree_lengths(master_path)
master_by_edge = {e: master_full[canon(e, ALL)] for e in full_splits}
root_kids = [frozenset(l.taxon.label for l in c.leaf_nodes()) for c in tree.seed_node.child_nodes()]
for e in root_kids:
    master_by_edge[e] = master_full[canon(e, ALL)] / len(root_kids)

# ---------- RER per gene ----------
rows_x, rows_y, rows_gene, rows_split = [], [], [], []
for f in done:
    tp, gl = tree_lengths(WORK / f / "tree.nwk", keep=ALL)
    if len(tp) < MIN_SP_GENE: continue
    mp = {}
    for e, L in master_by_edge.items():
        r = e & tp
        if 0 < len(r) < len(tp):
            k = canon(r, tp); mp[k] = mp.get(k, 0.0) + L
    keys = [k for k in gl if k in mp]
    rows_x += [mp[k] for k in keys]; rows_y += [gl[k] for k in keys]
    rows_gene += [f] * len(keys); rows_split += keys
df = pd.DataFrame({"gene": rows_gene, "x": np.sqrt(rows_x), "y": np.sqrt(rows_y)}); df["split"] = rows_split
res = np.empty(len(df))
for g, sub in df.groupby("gene").indices.items():
    xx, yy = df.x.values[sub], df.y.values[sub]
    A = np.column_stack([np.ones(len(xx)), xx]); res[sub] = yy - A @ np.linalg.lstsq(A, yy, rcond=None)[0]
df["res"] = res
df["bin"] = pd.qcut(df.x, 20, labels=False, duplicates="drop")
df["rer"] = df.res / df.groupby("bin").res.transform("std")
term = df[df.split.map(len) == 1].copy(); term["sp"] = term.split.map(lambda s: next(iter(s)))
R = term.pivot_table(index="gene", columns="sp", values="rer").reindex(columns=species)
R = R - R.median(axis=0)
R.to_csv(RES / "terminal_rer.tsv.gz", sep="\t", compression="gzip")
Rv = R.values; G = Rv.shape[0]
print(f"RER matrix {R.shape}, {np.isnan(Rv).mean():.1%} missing", flush=True)

# ---------- permulation-calibrated gene tests ----------
Ntips = len(tree.leaf_nodes())
height = {nd: (0.0 if nd.is_leaf() else (len(nd.leaf_nodes()) - 1) / (Ntips - 1)) for nd in tree.postorder_node_iter()}
C = np.eye(n)
for nd in tree.internal_nodes():
    kids = [[idx[l.taxon.label] for l in c.leaf_nodes()] for c in nd.child_nodes()]
    for a in range(len(kids)):
        for b in range(a + 1, len(kids)):
            C[np.ix_(kids[a], kids[b])] = 1 - height[nd]; C[np.ix_(kids[b], kids[a])] = 1 - height[nd]
Lc = np.linalg.cholesky(C + 1e-9 * np.eye(n))
mask = ~np.isnan(Rv)
ranks = np.full_like(Rv, np.nan)
for i in range(G): ranks[i, mask[i]] = stats.rankdata(Rv[i, mask[i]])


def auc_all(lab, sub=None):
    m = mask if sub is None else mask & sub[None, :]
    pos = m & (lab[None, :] == 1); n1 = pos.sum(1); n0 = m.sum(1) - n1
    rk = ranks if sub is None else np.where(m, ranks, np.nan)
    if sub is not None:  # re-rank within subset
        rk = np.full_like(Rv, np.nan)
        for i in range(G): rk[i, m[i]] = stats.rankdata(Rv[i, m[i]])
    rs = np.where(pos, rk, 0).sum(1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return (rs - n1 * (n1 + 1) / 2) / (n1 * n0)


def spear_all(yv):
    return np.array([stats.spearmanr(Rv[i, mask[i]], yv[mask[i]])[0] for i in range(G)])


def sim_trait():
    z = Lc @ rng.standard_normal(n); o = np.argsort(-z)
    ys = np.zeros(n); ys[o[:(y >= 1).sum()]] = 1; ys[o[:(y >= 2).sum()]] = 2; ys[o[:(y >= 3).sum()]] = 3
    return ys


def calib(obs, null):
    mu, sd = np.nanmean(null, 0), np.nanstd(null, 0); z = (obs - mu) / sd
    return z, 2 * stats.norm.sf(np.abs(z))


sims = [sim_trait() for _ in range(N_NULL)]
null_auc = np.array([auc_all((s >= 2).astype(int)) for s in sims])
null_rho = np.array([spear_all(s) for s in sims[:200]])
obs_auc, obs_rho = auc_all(ybin), spear_all(y)
z_auc, p_auc = calib(obs_auc, null_auc); z_rho, p_rho = calib(obs_rho, null_rho)
q = lambda p: multipletests(np.nan_to_num(p, nan=1), method="fdr_bh")[1]
out = pd.DataFrame({"family": R.index, "n_species": mask.sum(1), "auc_fractal": obs_auc, "z_auc": z_auc, "p_auc": p_auc,
                    "q_auc": q(p_auc), "spearman_order": obs_rho, "z_rho": z_rho, "p_rho": p_rho, "q_rho": q(p_rho)})


def subset_test(sub, label):
    ys_null = np.array([auc_all((s >= 2).astype(int), sub) for s in sims[:200]])
    o = auc_all(ybin, sub); z, p = calib(o, ys_null)
    out[f"z_auc_{label}"] = z; out[f"p_auc_{label}"] = p; out[f"q_auc_{label}"] = q(p)


subset_test(dtype == "genome", "genomes_only")
subset_test(~lowconf, "no_lowconf")
ga = pd.read_csv(PROC / "gene_family_assignments.tsv.gz", sep="\t")
ath = ga[(ga.sp == "ath") & ga.gf.isin(R.index)].groupby("gf").gene.apply(lambda s: ",".join(sorted(s)[:4]))
out["ath_genes"] = out.family.map(ath).fillna("")
out = out.sort_values("p_auc"); out.to_csv(RES / "rer_association.tsv", sep="\t", index=False)
null_hits = []
for k in range(min(20, N_NULL)):
    rest = np.delete(null_auc, k, 0)
    zz = (null_auc[k] - np.nanmean(rest, 0)) / np.nanstd(rest, 0)
    null_hits.append((q(2 * stats.norm.sf(np.abs(zz))) < 0.1).sum())

# ---------- leave-one-order-out classifier ----------
Xall = np.nan_to_num(Rv.T)


def cv(yb, groups=order):
    score = np.zeros(n)
    for c in np.unique(groups):
        te, tr = groups == c, groups != c
        if yb[tr].sum() == 0 or yb[tr].sum() == tr.sum(): continue
        a = np.empty(G)
        for i in range(G):
            m = mask[i] & tr; r = stats.rankdata(Rv[i, m]); l = yb[m]; n1, n0 = l.sum(), (1 - l).sum()
            a[i] = (r[l == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0) if n1 and n0 else 0.5
        top = np.argsort(-np.abs(a - 0.5))[:K]
        clf = LogisticRegression(C=0.1, class_weight="balanced", max_iter=2000).fit(Xall[tr][:, top], yb[tr])
        score[te] = clf.predict_proba(Xall[te][:, top])[:, 1]
    return score


score = cv(ybin)
auc, ap = roc_auc_score(ybin, score), average_precision_score(ybin, score)
null_cv = np.array([roc_auc_score((s >= 2).astype(int), cv((s >= 2).astype(int))) for s in sims[:N_NULL_CV]])
isgen = (dtype == "genome").astype(int)
auc_dtype = roc_auc_score(isgen, cv(isgen))
pred = pd.DataFrame({"sp": species, "species": st.loc[species, "species"].values, "order": order, "data_type": dtype,
                     "dissection": y, "fractal": ybin, "score": score}).sort_values("score", ascending=False)
pred["rank"] = np.arange(1, n + 1); pred.to_csv(RES / "rer_cv_predictions.tsv", sep="\t", index=False)

lines = [
    f"gene trees used: {G}; species {n} ({(dtype == 'genome').sum()} genomes, {(dtype != 'genome').sum()} transcriptomes); "
    f"fractal (order>=2) {ybin.sum()} in {len(set(order[ybin == 1]))} orders; nulls {N_NULL} (assoc), {N_NULL_CV} (cv)",
    "",
    "== Genes whose rate differs in fractal lineages (terminal-branch RER, permulation-calibrated) ==",
    f"q<0.1: {(out.q_auc < 0.1).sum()}  (random traits through the same pipeline: median {np.median(null_hits):.0f}, max {max(null_hits)})",
    f"q<0.1 correlation with dissection order: {(out.q_rho < 0.1).sum()}",
    f"q<0.1 genomes only: {(out.q_auc_genomes_only < 0.1).sum()}; without low-confidence scores: {(out.q_auc_no_lowconf < 0.1).sum()}",
    f"p<0.01: {(out.p_auc < 0.01).sum()} (expected {0.01 * G:.0f})",
    out.head(20)[["family", "n_species", "auc_fractal", "z_auc", "p_auc", "q_auc", "p_auc_genomes_only", "spearman_order", "ath_genes"]].to_string(index=False),
    "",
    "== Leave-one-order-out prediction from terminal RER ==",
    f"AUROC {auc:.3f}  avg precision {ap:.3f} (chance {ybin.mean():.3f})",
    f"null (Brownian-threshold traits): mean {null_cv.mean():.3f}, 95th pct {np.quantile(null_cv, .95):.3f}; P(null >= obs) = {(null_cv >= auc).mean():.3f}",
    f"artifact check: AUROC for predicting genome vs transcriptome from the same features = {auc_dtype:.3f}",
    pred[pred.fractal == 1].head(40).to_string(index=False),
]
(RES / "summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
