"""Per-family protein alignments and branch lengths on the fixed species topology.

Families: backbone families present in >= 80% of included species and >= 7/9
fractal species, median <= 2 copies. For each family:
  1. MAFFT-align all member proteins from included species.
  2. Keep one gene per species: the copy with highest identity to the column
     consensus (reduces paralog/homeolog noise).
  3. trimAl -gappyout.
  4. IQ-TREE (LG+G4) optimizes branch lengths on the species tree pruned to
     the species present (-te: topology fixed).
Usage: python3 06_gene_trees.py [n_families|all] [n_workers]
"""
import gzip, pathlib, subprocess, sys, collections
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd, dendropy

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC = ROOT / "data/processed"
WORK = ROOT / "data/work/gene_trees"; WORK.mkdir(parents=True, exist_ok=True)
BIN = pathlib.Path("/home/user/tools/env/bin")
LIMIT = sys.argv[1] if len(sys.argv) > 1 else "all"
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 4


def select_families(inc, fractal):
    m = pd.read_csv(PROC / "family_copy_number.tsv.gz", sep="\t", index_col=0).drop(columns="_total_genes").loc[inc]
    pres = (m > 0).mean(); med = m.where(m > 0).median()
    pf = (m.loc[fractal] > 0).sum()
    return m.columns[(pres >= 0.8) & (med <= 2) & (pf >= 7)].tolist()


def load_proteins(genes_by_sp):
    seqs = {}
    for sp, genes in genes_by_sp.items():
        f = next((ROOT / "data/raw/proteomes").glob(f"{sp}.*gz"))
        want, cur, buf = set(genes), None, []
        with gzip.open(f, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    if cur in want: seqs[(sp, cur)] = "".join(buf).rstrip("*")
                    h = line[1:].strip()
                    cur = h.split("|")[1].strip() if "|" in h else h.split()[0]
                    buf = []
                else:
                    buf.append(line.strip())
            if cur in want: seqs[(sp, cur)] = "".join(buf).rstrip("*")
    return seqs


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        if line.startswith(">"): name = line[1:].strip(); out[name] = []
        else: out[name].append(line.strip())
    return {k: "".join(v) for k, v in out.items()}


def run_family(fam):
    d = WORK / fam
    if (d / "done").exists():
        return fam, "cached"
    try:
        aln = read_fasta(d / "all.aln") if (d / "all.aln").exists() else None
        if aln is None:
            with open(d / "all.aln", "w") as out:
                subprocess.run([BIN / "mafft", "--retree", "2", "--thread", "1", "--quiet", d / "members.fa"],
                               stdout=out, check=True)
            aln = read_fasta(d / "all.aln")
        names = list(aln); A = np.array([list(aln[n].upper()) for n in names])
        cons = []
        for col in A.T:
            v, c = np.unique(col[col != "-"], return_counts=True)
            cons.append(v[c.argmax()] if len(v) else "-")
        cons = np.array(cons)
        ident = ((A == cons) & (A != "-")).sum(1) / max((cons != "-").sum(), 1)
        best = {}
        for n, s in zip(names, ident):
            sp = n.split("__")[0]
            if sp not in best or s > best[sp][1]: best[sp] = (n, s)
        keep = [best[sp][0] for sp in best]
        B = A[[names.index(k) for k in keep]]
        B = B[:, (B != "-").any(0)]
        with open(d / "rep.aln", "w") as fh:
            for k, row in zip(keep, B): fh.write(f">{k.split('__')[0]}\n{''.join(row)}\n")
        subprocess.run([BIN / "trimal", "-in", d / "rep.aln", "-out", d / "rep.trim.aln", "-gappyout"],
                       check=True, capture_output=True)
        trim = {k: v for k, v in read_fasta(d / "rep.trim.aln").items() if v.strip("-X")}
        if len(trim) < 20:
            raise ValueError(f"only {len(trim)} sequences after trimming")
        with open(d / "rep.trim.aln", "w") as fh:  # trimAl can leave all-gap rows or drop rows
            for k, v in trim.items(): fh.write(f">{k}\n{v}\n")
        t = dendropy.Tree.get(path=PROC / "species_tree.nwk", schema="newick", preserve_underscores=True)
        t.retain_taxa_with_labels(list(trim)); t.suppress_unifurcations()
        (d / "topology.nwk").write_text(t.as_string(schema="newick").replace("[&R] ", ""))
        subprocess.run([BIN / "iqtree", "-s", d / "rep.trim.aln", "-te", d / "topology.nwk", "-m", "LG+G4",
                        "-nt", "1", "-quiet", "-redo", "-pre", d / "iq"], check=True, capture_output=True)
        (d / "done").write_text("ok")
        for f in ["iq.log", "iq.ckp.gz", "iq.uniqueseq.phy", "iq.model.gz", "all.aln", "members.fa"]:
            (d / f).unlink(missing_ok=True)
        return fam, "ok"
    except Exception as e:
        return fam, f"error {e}"


if __name__ == "__main__":
    ph = pd.read_csv(ROOT / "data/phenotypes/leaf_dissection_curated.tsv", sep="\t", comment="#")
    inc = ph[ph.include == 1].sp.tolist()
    fractal = ph[(ph.include == 1) & (ph.dissection_order >= 2)].sp.tolist()
    fams = select_families(inc, fractal)
    rng = np.random.default_rng(0)
    fams = list(rng.permutation(fams))  # random order: partial runs are an unbiased subset
    if LIMIT != "all":
        fams = fams[:int(LIMIT)]
    todo = [f for f in fams if not (WORK / f / "done").exists()]
    print(f"{len(fams)} families, {len(todo)} to run", flush=True)
    ga = pd.read_csv(PROC / "gene_family_assignments.tsv.gz", sep="\t")
    ga = ga[ga.gf.isin(todo) & ga.sp.isin(inc)]
    seqs = load_proteins(ga.groupby("sp").gene.apply(list).to_dict())
    for fam, sub in ga.groupby("gf"):
        d = WORK / fam; d.mkdir(exist_ok=True)
        if not (d / "members.fa").exists():
            with open(d / "members.fa", "w") as fh:
                for sp, g in zip(sub.sp, sub.gene):
                    s = seqs.get((sp, g))
                    if s: fh.write(f">{sp}__{g}\n{s}\n")
    (PROC / "rer_families.txt").write_text("\n".join(fams) + "\n")
    status = collections.Counter()
    with ProcessPoolExecutor(NW) as ex:
        for i, (fam, st) in enumerate(ex.map(run_family, todo, chunksize=1)):
            status[st.split()[0]] += 1
            if st.startswith("error"): print(fam, st[:300], flush=True)
            if i % 100 == 0: print(i, dict(status), flush=True)
    print("finished", dict(status))
