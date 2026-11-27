"""Branch-level rate tests with ancestral reconstruction and a trait-matched null.

Fixes two problems of 14_rer_expanded.py:
  * only terminal branches were tested: here every branch of every gene tree is used
    (RERconverge-style), with the trait reconstructed on internal branches, plus a
    separate test of the branches where the trait arose (transitions);
  * the Brownian-threshold null was more clumped than the real trait: here null traits
    are simulated under a 2-state Mk model fitted to the real trait on the master tree,
    and kept only if their number of trait tips and of reconstructed origins match the
    observed ones (+-25% and +-50%).

Trait: a 0/1 column of data/phenotypes/traits.tsv (NA = excluded from that trait).
RER per gene branch (sqrt-length regression on the master tree, binned scaling, terminal
branches centered per species) is trait independent and cached in data/work/branch_rer.npz.
Per gene: AUC of RER on foreground (trait-state) branches vs background, and on
transition branches vs background; calibrated against the matched null.
Prediction: terminal RER, leave-one-origin-out (each reconstructed origin is a fold).
Outputs: results/branch_tests/<trait>/{summary.txt, genes.tsv, cv_predictions.tsv}
Usage: python3 16_branch_tests.py <trait> [n_null] [K]
"""
import pathlib, sys, re, os
import numpy as np, pandas as pd
from scipy import optimize, sparse, stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC, RESX = ROOT / "data/processed", ROOT / "results/rer_expanded"
# BT_GENE_DIR / BT_TAG: run the same tests on another set of gene trees (e.g. regulator subclades)
TAG = os.environ.get("BT_TAG", "")
WORK = pathlib.Path(os.environ.get("BT_GENE_DIR", ROOT / "data/work/gene_trees_x"))
CACHE = ROOT / f"data/work/branch_rer{TAG}.npz"
TRAIT = sys.argv[1] if len(sys.argv) > 1 else "fractal_old"
N_NULL = int(sys.argv[2]) if len(sys.argv) > 2 else 500
K = int(sys.argv[3]) if len(sys.argv) > 3 else 50
MIN_SP_GENE = int(os.environ.get("BT_MIN_SP", 100))
OUT = ROOT / "results/branch_tests" / (TRAIT + TAG); OUT.mkdir(parents=True, exist_ok=True)
rng = np.random.default_rng(3)


# ---------------- minimal newick handling with tip bitmasks ----------------
def parse_newick(s):
    """Return list of nodes: dict(parent, children, length, label)."""
    s = s.strip().rstrip(";")
    nodes, stack, i, cur = [], [], 0, None
    tok = re.compile(r"\(|\)|,|:[^,();]+|[^,():;]+")
    new = lambda p: (nodes.append({"parent": p, "children": [], "length": None, "label": None}), len(nodes) - 1)[1]
    root = new(-1); cur = root
    for m in tok.finditer(s):
        t = m.group()
        if t == "(":
            c = new(cur); nodes[cur]["children"].append(c); stack.append(cur); cur = c
        elif t == ",":
            p = stack[-1]; c = new(p); nodes[p]["children"].append(c); cur = c
        elif t == ")":
            cur = stack.pop()
        elif t.startswith(":"):
            nodes[cur]["length"] = float(t[1:])
        else:
            nodes[cur]["label"] = t
    return nodes


def masks(nodes, tip_index):
    """Bitmask of tips (indices in tip_index) below each node (postorder accumulation)."""
    m = [0] * len(nodes)
    for i in range(len(nodes) - 1, -1, -1):  # children always have larger index than parent
        n = nodes[i]
        if not n["children"]:
            lab = n["label"]
            m[i] = (1 << tip_index[lab]) if lab in tip_index else 0
        for c in n["children"]:
            m[i] |= m[c]
    return m


def canon(r, tp):
    low = tp & -tp
    return tp ^ r if r & low else r


# ---------------- species, full rooted tree, master lengths ----------------
sp_list = pd.read_csv(RESX / "rer_cv_predictions.tsv", sep="\t").sp.tolist()  # the 1,133 species of step 14
species = sorted(sp_list); n = len(species); tix = {s: i for i, s in enumerate(species)}
FULL = (1 << n) - 1
topo = parse_newick((RESX / "master_topology.nwk").read_text())
tmask = masks(topo, tix)
master_nodes = parse_newick((RESX / "master.nwk").read_text())
mm = masks(master_nodes, tix)
master_len = {canon(mm[i], FULL): (nd["length"] or 0.0) for i, nd in enumerate(master_nodes) if nd["parent"] >= 0}
edges = [i for i, nd in enumerate(topo) if nd["parent"] >= 0]           # full-tree edges = non-root nodes
emask = np.array([tmask[i] for i in edges], dtype=object)
root_kids = topo[0]["children"]
elen = np.array([master_len[canon(tmask[i], FULL)] / (len(root_kids) if topo[i]["parent"] == 0 else 1) for i in edges])
E = len(edges); edge_of_node = {nd: k for k, nd in enumerate(edges)}

# ---------------- RER per gene branch (cached) ----------------
if CACHE.exists():
    z = np.load(CACHE, allow_pickle=True)
    gene_names = z["gene_names"]; row_gene, rer, row_tip = z["row_gene"], z["rer"], z["row_tip"]
    W = sparse.csr_matrix((z["w_data"], z["w_ind"], z["w_ptr"]), shape=tuple(z["w_shape"]))
else:
    fams = sorted(p.parent.name for p in WORK.glob("*/done"))
    gene_names, rg, xs, ys, tipi, wr, wc, wv = [], [], [], [], [], [], [], []
    row = 0
    for f in fams:
        gn = parse_newick((WORK / f / "tree.nwk").read_text())
        gmask = masks(gn, tix); tp = gmask[0]
        if bin(tp).count("1") < MIN_SP_GENE: continue
        gl = {}
        for i, nd in enumerate(gn):
            if nd["parent"] < 0: continue
            k = canon(gmask[i], tp)
            if 0 < k < tp and k != tp: gl[k] = gl.get(k, 0.0) + (nd["length"] or 0.0)
        mp = {}
        for e in range(E):
            r = emask[e] & tp
            if r and r != tp:
                mp.setdefault(canon(r, tp), []).append(e)
        gi = len(gene_names); gene_names.append(f)
        for k, L in gl.items():
            if k not in mp: continue
            es = mp[k]
            xs.append(elen[es].sum()); ys.append(L); rg.append(gi)
            pc = bin(k).count("1")
            if pc == 1: tipi.append(k.bit_length() - 1)
            elif bin(tp ^ k).count("1") == 1: tipi.append((tp ^ k).bit_length() - 1)
            else: tipi.append(-1)
            wr += [row] * len(es); wc += es; wv += list(elen[es]); row += 1
        if gi % 500 == 0: print("rer", gi, row, flush=True)
    row_gene, row_tip = np.array(rg), np.array(tipi)
    x, y = np.sqrt(xs), np.sqrt(ys)
    res = np.empty(len(x))
    starts = np.r_[0, np.flatnonzero(np.diff(row_gene)) + 1, len(x)]
    for a, b in zip(starts[:-1], starts[1:]):
        A = np.column_stack([np.ones(b - a), x[a:b]]); res[a:b] = y[a:b] - A @ np.linalg.lstsq(A, y[a:b], rcond=None)[0]
    bins = pd.qcut(x, 20, labels=False, duplicates="drop")
    sd = pd.Series(res).groupby(bins).transform("std").values
    rer = res / sd
    # center terminal branches per species (removes species-level offsets such as data type)
    term = row_tip >= 0
    med = pd.Series(rer[term]).groupby(row_tip[term]).median()
    rer[term] -= med.reindex(row_tip[term]).values
    gene_names = np.array(gene_names)
    W = sparse.csr_matrix((wv, (wr, wc)), shape=(row, E))
    np.savez(CACHE, gene_names=gene_names, row_gene=row_gene, rer=rer, row_tip=row_tip,
             w_data=W.data, w_ind=W.indices, w_ptr=W.indptr, w_shape=np.array(W.shape))
G = len(gene_names); R_ = len(rer)
wsum = np.asarray(W.sum(1)).ravel()
order_rows = np.lexsort((rer, row_gene))           # rows sorted by gene, then RER
seg = np.r_[0, np.flatnonzero(np.diff(row_gene[order_rows])) + 1]
print(f"{G} genes, {R_} gene branches, {E} species-tree edges", flush=True)

# ---------------- trait and Mk reconstruction on the pruned tree ----------------
traits = pd.read_csv(ROOT / "data/phenotypes/traits.tsv", sep="\t").set_index("sp")
tv = traits[TRAIT].reindex(species)
known = tv.notna().values; yt = np.where(known, tv.fillna(0).values, -1).astype(int)
S = sum(1 << i for i in range(n) if known[i])

# pruned tree over known tips: node = full-tree node whose mask&S is a distinct non-empty clade
clade_of = {}
for i in range(len(topo)):
    r = tmask[i] & S
    if r: clade_of.setdefault(r, []).append(i)   # several full nodes may share one pruned clade (a path)
pr_clades = sorted(clade_of, key=lambda r: -bin(r).count("1"))  # parents before children
pidx = {r: j for j, r in enumerate(pr_clades)}
P = len(pr_clades); par = np.full(P, -1); blen = np.zeros(P)
for j, r in enumerate(pr_clades):
    top = min(clade_of[r])                          # highest full node of the path
    nodes_on_path = clade_of[r]
    blen[j] = sum(elen[edge_of_node[i]] for i in nodes_on_path if i in edge_of_node)
    p = topo[top]["parent"]
    while p >= 0 and (tmask[p] & S) == r: p = topo[p]["parent"]
    if p >= 0: par[j] = pidx[tmask[p] & S]
blen = np.maximum(blen, 1e-6)
kids = [[] for _ in range(P)]
for j in range(P):
    if par[j] >= 0: kids[par[j]].append(j)
is_tip = np.array([not k for k in kids])
tip_sp = {j: (pr_clades[j].bit_length() - 1) for j in range(P) if is_tip[j]}
post = sorted(range(P), key=lambda j: bin(pr_clades[j]).count("1"))   # children before parents
pre = post[::-1]
root_p = int(np.flatnonzero(par < 0)[0])
# full-tree edge -> pruned node whose path contains it (-1: no known tips below / above all of them)
edge_pr = np.full(E, -1)
for j, r in enumerate(pr_clades):
    if par[j] < 0: continue
    for i in clade_of[r]:
        if i in edge_of_node: edge_pr[edge_of_node[i]] = j
has_pr = edge_pr >= 0
par_safe = np.where(par >= 0, par, 0)


def pmat(t, q01, q10):
    s = q01 + q10; e = np.exp(-s * t)[:, None, None]
    pi0, pi1 = q10 / s, q01 / s
    stat = np.array([[pi0, pi1], [pi0, pi1]])
    return stat[None] + e * np.array([[pi1, -pi1], [-pi0, pi0]])[None]


def prune(tips_state, q01, q10, Pm=None):
    Pm = pmat(blen, q01, q10) if Pm is None else Pm
    Lk = np.ones((P, 2)); logsc = 0.0
    for j in post:
        if is_tip[j]:
            st = tips_state[tip_sp[j]]; Lk[j] = [st == 0, st == 1]
        else:
            v = np.ones(2)
            for c in kids[j]: v *= Pm[c] @ Lk[c]
            m = v.max(); Lk[j] = v / m; logsc += np.log(m)
    pi = np.array([q10, q01]) / (q01 + q10)
    return Lk, logsc + np.log(pi @ Lk[root_p]), Pm, pi


def fit(tips_state):
    f = lambda lq: -prune(tips_state, *np.exp(lq))[1]
    r = optimize.minimize(f, np.log([0.5, 0.5]), method="Nelder-Mead", options={"xatol": 1e-3, "fatol": 1e-3})
    return np.exp(r.x)


def marginal(tips_state, q01, q10):
    """Posterior P(state=1) at every pruned node (up-down algorithm)."""
    Lk, _, Pm, pi = prune(tips_state, q01, q10)
    up = np.zeros((P, 2)); up[root_p] = pi
    post1 = np.zeros(P)
    for j in pre:  # preorder
        if j != root_p:
            p = par[j]
            sib = np.ones(2)
            for c in kids[p]:
                if c != j: sib *= Pm[c] @ Lk[c]
            a = up[p] * sib; a /= a.sum()
            up[j] = a @ Pm[j]; up[j] /= up[j].sum()
        b = up[j] * Lk[j]; post1[j] = b[1] / b.sum()
    return post1


def edge_states(post1):
    """Map pruned-node posteriors to full-tree edges: fg (state 1), transition (0->1 on the path), defined."""
    st = post1 > 0.5
    trans_node = st & ~st[par_safe] & (par >= 0)
    j = np.where(has_pr, edge_pr, 0)
    fg = np.where(has_pr, st[j], 0).astype(float); tr = np.where(has_pr, trans_node[j], 0).astype(float)
    return fg, tr, has_pr.astype(float), int(trans_node.sum())


def sim_tips(q01, q10, Pm, pi):
    st = np.zeros(P, int); st[root_p] = rng.random() < pi[1]
    u = rng.random(P)
    for j in pre:
        if j != root_p: st[j] = u[j] < Pm[j][st[par[j]], 1]
    out = np.full(n, -1)
    for j, s in tip_sp.items(): out[s] = st[j]
    return out


# ---------------- per-gene statistics ----------------
def gene_auc(pos, neg):
    """AUC per gene of RER for pos rows vs neg rows (bool arrays over rows)."""
    v = (pos | neg)[order_rows]; p = pos[order_rows]
    rk = np.cumsum(v); rk = rk - np.repeat(np.r_[0, rk][seg], np.diff(np.r_[seg, len(rk)]))
    n1 = np.add.reduceat(p.astype(float), seg); nv = np.add.reduceat(v.astype(float), seg); n0 = nv - n1
    rs = np.add.reduceat(np.where(p, rk, 0).astype(float), seg)
    with np.errstate(invalid="ignore", divide="ignore"):
        auc = (rs - n1 * (n1 + 1) / 2) / (n1 * n0)
        # Mann-Whitney z: removes the dependence of the AUC's spread on n1/n0, which differ between
        # null traits; remaining excess spread (branches are not independent) is calibrated by the null
        z = (auc - 0.5) / np.sqrt((n1 + n0 + 1) / (12 * n1 * n0))
    return auc, n1, z


def branch_stats(fg, tr, dfn):
    ok = np.isclose(W @ dfn, wsum)
    fgw = (W @ fg) / wsum; isfg = ok & (fgw > 0.5); isbg = ok & (fgw <= 0.5)
    istr = ok & ((W @ tr) > 0)
    a_fg, n_fg, z_fg = gene_auc(isfg, isbg)
    a_tr, n_tr, z_tr = gene_auc(istr, isbg & ~istr)
    return a_fg, a_tr, n_fg, n_tr, z_fg, z_tr


q01, q10 = fit(yt)
post_obs = marginal(yt, q01, q10)
fg, tr, dfn, n_orig = edge_states(post_obs)
a_fg, a_tr, n_fg, n_tr, mw_fg, mw_tr = branch_stats(fg, tr, dfn)
n1_obs = int((yt == 1).sum())
print(f"trait {TRAIT}: {n1_obs} of {known.sum()} species; Mk q01={q01:.3f} q10={q10:.3f}; {n_orig} reconstructed origins", flush=True)

# matched null
Pm = pmat(blen, q01, q10); pi = np.array([q10, q01]) / (q01 + q10)
null_fg, null_tr, null_tips, tries = [], [], [], 0
while len(null_fg) < N_NULL and tries < 200 * N_NULL:
    tries += 1
    ys_ = sim_tips(q01, q10, Pm, pi); k1 = (ys_ == 1).sum()
    if not (0.75 * n1_obs <= k1 <= 1.25 * n1_obs): continue
    ps = marginal(ys_, q01, q10); f_, t_, d_, no = edge_states(ps)
    if not (0.5 * n_orig <= no <= 1.5 * n_orig): continue
    _, _, _, _, a, b = branch_stats(f_, t_, d_); null_fg.append(a); null_tr.append(b); null_tips.append(ys_)
    if len(null_fg) % 50 == 0: print("null", len(null_fg), "tries", tries, flush=True)
null_fg, null_tr = np.array(null_fg), np.array(null_tr)
np.savez(OUT / "null_traits.npz", species=np.array(species), observed=yt, null=np.array(null_tips))  # reused by 17_regulators.py


def calib(obs, null):
    mu, sd = np.nanmean(null, 0), np.nanstd(null, 0); z = (obs - mu) / sd
    return z, 2 * stats.norm.sf(np.abs(z))


q = lambda p: multipletests(np.nan_to_num(p, nan=1), method="fdr_bh")[1]
z_fg, p_fg = calib(mw_fg, null_fg); z_tr, p_tr = calib(mw_tr, null_tr)  # nulls hold Mann-Whitney z
null_hits = []
for k in range(min(20, len(null_fg))):
    rest = np.delete(null_fg, k, 0); zz = (null_fg[k] - np.nanmean(rest, 0)) / np.nanstd(rest, 0)
    null_hits.append((q(2 * stats.norm.sf(np.abs(zz))) < 0.1).sum())
ga = pd.read_csv(PROC / "gene_family_assignments.tsv.gz", sep="\t")
ath = ga[(ga.sp == "ath") & ga.gf.isin(gene_names)].groupby("gf").gene.apply(lambda s: ",".join(sorted(s)[:4]))
out = pd.DataFrame({"family": gene_names, "n_fg_branches": n_fg, "auc_fg": a_fg, "z_fg": z_fg, "p_fg": p_fg, "q_fg": q(p_fg),
                    "n_transition_branches": n_tr, "auc_transition": a_tr, "z_tr": z_tr, "p_tr": p_tr, "q_tr": q(p_tr)})
out["ath_genes"] = out.family.map(ath).fillna("")
out = out.sort_values("p_fg"); out.to_csv(OUT / "genes.tsv", sep="\t", index=False)

# ---------------- leave-one-origin-out prediction from terminal RER ----------------
term = row_tip >= 0
Rm = pd.DataFrame({"g": row_gene[term], "s": row_tip[term], "r": rer[term]}).pivot_table(index="g", columns="s", values="r")
Rm = Rm.reindex(index=range(G), columns=range(n)); Rv = Rm.values; mask = ~np.isnan(Rv)
X = np.nan_to_num(Rv.T)


def folds(post1, ytips):
    """Fold per species: the reconstructed origin clade containing it; species outside any origin: by NCBI order."""
    st = post1 > 0.5; f = np.array(["o:" + str(o) for o in orders], dtype=object)
    for j in pre:
        if par[j] >= 0 and st[j] and not st[par[j]]:
            r = pr_clades[j]
            for i in range(n):
                if r >> i & 1: f[i] = f"origin{j}"
    return f


st_tab = pd.read_csv(PROC / "species_table.tsv", sep="\t").set_index("sp")
orders = st_tab.loc[species, "order"].fillna("NA").values


ORD_T = np.argsort(np.where(mask, Rv, np.inf), axis=1); MASK_T = np.take_along_axis(mask, ORD_T, 1)


def tip_auc(lab, sub):
    ms = MASK_T & sub[ORD_T]; rk = np.cumsum(ms, 1); pos = ms & (lab[ORD_T] == 1)
    n1 = pos.sum(1); n0 = ms.sum(1) - n1
    with np.errstate(invalid="ignore", divide="ignore"):
        a = (np.where(pos, rk, 0).sum(1) - n1 * (n1 + 1) / 2) / (n1 * n0)
    return np.nan_to_num(a, nan=0.5)


def cv(lab, fo):
    score = np.full(n, np.nan); kn = lab >= 0
    for c in np.unique(fo[kn]):
        te = kn & (fo == c); trn = kn & (fo != c)
        if lab[trn].sum() in (0, trn.sum()): continue
        a = tip_auc(lab, trn); top = np.argsort(-np.abs(a - 0.5))[:K]
        clf = LogisticRegression(C=0.1, class_weight="balanced", max_iter=2000).fit(X[trn][:, top], lab[trn])
        score[te] = clf.predict_proba(X[te][:, top])[:, 1]
    m = kn & ~np.isnan(score)
    return roc_auc_score(lab[m], score[m]), score


auc_cv, score = cv(yt, folds(post_obs, yt))
null_cv = []
for ys_ in null_tips[:min(30, len(null_tips))]:
    null_cv.append(cv(ys_, folds(marginal(ys_, q01, q10), ys_))[0])
null_cv = np.array(null_cv)
pred = pd.DataFrame({"sp": species, "species": st_tab.loc[species, "species"].values, "order": orders,
                     "trait": yt, "score": score}).sort_values("score", ascending=False)
pred.to_csv(OUT / "cv_predictions.tsv", sep="\t", index=False)

lines = [
    f"trait {TRAIT}: {n1_obs} of {known.sum()} species with state 1; {n_orig} reconstructed origins (Mk q01={q01:.3f}, q10={q10:.3f})",
    f"genes {G}; gene branches {R_}; matched nulls {len(null_fg)} (of {tries} simulations)",
    "",
    "== Rate on trait-state branches (all branches, terminal + internal) ==",
    f"q<0.1: {(out.q_fg < 0.1).sum()} (null traits: median {np.median(null_hits):.0f}, max {max(null_hits) if null_hits else 'NA'}); "
    f"p<0.01: {(out.p_fg < 0.01).sum()} (expected {0.01 * G:.0f}); p<0.001: {(out.p_fg < 0.001).sum()} (expected {0.001 * G:.0f})",
    "== Rate on branches where the trait arose (transitions) ==",
    f"q<0.1: {(out.q_tr < 0.1).sum()}; p<0.01: {(out.p_tr < 0.01).sum()} (expected {0.01 * G:.0f}); p<0.001: {(out.p_tr < 0.001).sum()}",
    "",
    out.head(15).round(4).to_string(index=False),
    "",
    out.sort_values("p_tr").head(10).round(4).to_string(index=False),
    "",
    "== Leave-one-origin-out prediction (terminal RER) ==",
    f"AUROC {auc_cv:.3f}; matched-null traits: mean {null_cv.mean():.3f}, 95th pct {np.quantile(null_cv, .95):.3f}, "
    f"P(null >= obs) = {(null_cv >= auc_cv).mean():.3f} (n={len(null_cv)})",
]
(OUT / "summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines))
