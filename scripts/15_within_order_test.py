"""Within-order check of the expanded RER results.

Compares fractal vs non-fractal species only within the same NCBI order
(orders containing both), so between-order rate differences (e.g. ferns vs
seed plants) cannot create associations. Statistic per gene: sum over fractal
species of the within-order centered rank of terminal RER. Null: fractal labels
shuffled within each order (counts per order preserved), N_PERM times.
Output: results/rer_expanded/within_order.tsv, within_order_summary.txt
Usage: python3 15_within_order_test.py [n_perm]
"""
import pathlib, sys
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results/rer_expanded"
N_PERM = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
rng = np.random.default_rng(11)

R = pd.read_csv(RES / "terminal_rer.tsv.gz", sep="\t", index_col=0)
pr = pd.read_csv(RES / "rer_cv_predictions.tsv", sep="\t").set_index("sp").loc[R.columns]
y, order, dtype = pr.fractal.values, pr.order.fillna("NA").values, pr.data_type.values
mixed = [o for o in np.unique(order) if o != "NA" and 0 < y[order == o].sum() < (order == o).sum()]
use = np.isin(order, mixed)


def run(sub, label):
    keep = use & sub
    Rv, yy, oo = R.values[:, keep], y[keep], order[keep]
    U = np.zeros_like(Rv)
    for o in np.unique(oo):
        j = np.where(oo == o)[0]
        for i in range(Rv.shape[0]):
            m = ~np.isnan(Rv[i, j])
            if m.sum() > 1 and 0 < yy[j][m].sum() < m.sum():
                U[i, j[m]] = (stats.rankdata(Rv[i, j[m]]) - 0.5) / m.sum() - 0.5
    T = U @ yy
    P = np.zeros((len(yy), N_PERM))
    for o in np.unique(oo):
        j = np.where(oo == o)[0]
        for k in range(N_PERM): P[j, k] = rng.permutation(yy[j])
    Tn = U @ P
    z = (T - Tn.mean(1)) / Tn.std(1)
    p = (1 + (np.abs(Tn - Tn.mean(1, keepdims=True)) >= np.abs(T - Tn.mean(1))[:, None]).sum(1)) / (N_PERM + 1)
    q = multipletests(p, method="fdr_bh")[1]
    # same pipeline on 20 permuted label sets: how many q<0.1 by chance
    null_hits = []
    for k in range(20):
        zk = (Tn[:, k] - Tn.mean(1)) / Tn.std(1)
        null_hits.append((multipletests(2 * stats.norm.sf(np.abs(zk)), method="fdr_bh")[1] < 0.1).sum())
    n_ord = len(np.unique(oo[yy == 1]))
    return pd.DataFrame({f"z_{label}": z, f"p_{label}": p, f"q_{label}": q}, index=R.index), \
        f"{label}: {keep.sum()} species, {yy.sum()} fractal in {n_ord} orders; q<0.1: {(q < 0.1).sum()} " \
        f"(permuted labels: median {np.median(null_hits):.0f}, max {max(null_hits)}); p<0.01: {(p < 0.01).sum()} " \
        f"(expected {0.01 * len(p):.0f}); p<0.001: {(p < 0.001).sum()} (expected {0.001 * len(p):.0f})"


out, lines = [], []
for sub, label in [(np.ones(len(y), bool), "all"), (dtype == "genome", "genomes_only"), (dtype != "genome", "transcriptomes_only")]:
    d, line = run(sub, label); out.append(d); lines.append(line); print(line, flush=True)
df = pd.concat(out, axis=1).sort_values("p_all")
assoc = pd.read_csv(RES / "rer_association.tsv", sep="\t").set_index("family")
df["ath_genes"] = assoc.ath_genes.reindex(df.index)
df.to_csv(RES / "within_order.tsv", sep="\t")
top = df.head(15)
lines += ["", f"orders with both fractal and non-fractal species: {len(mixed)}", "", top.round(4).to_string()]
(RES / "within_order_summary.txt").write_text("\n".join(lines) + "\n")
print("\n".join(lines[3:]))
