"""Robustness checks for regulator copy-number associations (results/regulators/).

For each subclade x trait: Mann-Whitney z of log2(copies+1) relative to the species'
background duplication level (mean log2(copies+1) over the other regulator subclades),
calibrated against the matched null traits of 16_branch_tests.py, in
  genomes / 1KP transcriptomes (independent species) / genomes without the leading families;
plus a within-order permutation test (genomes, orders containing both states).
Output: results/regulators/checks.tsv and checks_summary.txt
Usage: python3 18_regulator_checks.py [n_perm]
"""
import pathlib, sys
import numpy as np, pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

ROOT = pathlib.Path(__file__).resolve().parents[1]
RES = ROOT / "results/regulators"
N_PERM = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
TRAITS = ["dissected_angio", "compound_angio", "divided_angio", "fractal_angio", "fractal"]
cn = pd.read_csv(RES / "copy_number.tsv", sep="\t").set_index("sp")
st = pd.read_csv(ROOT / "data/processed/species_table.tsv", sep="\t").set_index("sp")
rng = np.random.default_rng(5)
rows = []
for trait in TRAITS:
    z = np.load(ROOT / f"results/branch_tests/{trait}/null_traits.npz")
    species = list(z["species"]); obs, null = z["observed"], z["null"]
    c = cn.reindex(species).fillna(0); L = np.log2(c + 1)
    dtype = st.loc[species, "data_type"].values; fam = st.loc[species, "family"].values
    order = st.loc[species, "order"].fillna("NA").values
    genome = dtype == "genome"
    for sub in c.columns:
        x = (L[sub] - L.drop(columns=[sub]).mean(axis=1)).values

        def mwz(lab, m0):
            m = m0 & (lab >= 0); y = lab[m]; n1 = y.sum(); n0 = len(y) - n1
            if n1 < 3 or n0 < 3: return np.nan
            r = stats.rankdata(x[m])
            return ((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0) - 0.5) / np.sqrt((n1 + n0 + 1) / (12 * n1 * n0))

        def calib(m0):
            zo = mwz(obs, m0); zn = np.array([mwz(v, m0) for v in null])
            zc = (zo - np.nanmean(zn)) / np.nanstd(zn); return zc, 2 * stats.norm.sf(abs(zc))
        zg, pg = calib(genome); zt, pt = calib(dtype == "transcriptome")
        # leave out the 3 families contributing most trait-state genomes with above-median x
        m = genome & (obs == 1)
        top = pd.Series(fam[m][x[m] > np.median(x[genome & (obs >= 0)])]).value_counts().index[:3].tolist()
        zl, pl = calib(genome & ~np.isin(fam, top))
        mg = genome & (obs >= 0)
        mixed = [o for o in np.unique(order[mg]) if 0 < obs[mg & (order == o)].sum() < (mg & (order == o)).sum()]
        idx = [np.flatnonzero(mg & (order == o)) for o in mixed]
        rk = [(stats.rankdata(x[j]) - 0.5) / len(j) - 0.5 for j in idx]
        T = sum(r[obs[j] == 1].sum() for r, j in zip(rk, idx))
        Tn = np.array([sum(r[rng.permutation(obs[j]) == 1].sum() for r, j in zip(rk, idx)) for _ in range(N_PERM)])
        pw = (1 + (np.abs(Tn) >= abs(T)).sum()) / (N_PERM + 1)
        rows.append((trait, sub, zg, pg, zt, pt, ";".join(top), zl, pl, len(mixed), pw))
res = pd.DataFrame(rows, columns=["trait", "subclade", "z_genomes", "p_genomes", "z_transcriptomes", "p_transcriptomes",
                                  "left_out_families", "z_genomes_leftout", "p_genomes_leftout", "n_orders", "p_within_order"])
res["q_genomes"] = multipletests(res.p_genomes.fillna(1), method="fdr_bh")[1]  # across all traits x subclades
res = res.sort_values("p_genomes"); res.to_csv(RES / "checks.tsv", sep="\t", index=False)
lines = [f"{len(res)} tests (5 traits x {cn.shape[1]} subclades); copy number relative to background duplication level",
         f"q<0.1 (genomes, BH over all tests): {(res.q_genomes < 0.1).sum()}", "",
         res.head(12).round(4).to_string(index=False)]
(RES / "checks_summary.txt").write_text("\n".join(lines) + "\n"); print("\n".join(lines))
