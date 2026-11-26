"""Expanded gene trees: add NCBI/1KP candidates to the PLAZA alignments, then
branch lengths on the fixed expanded species topology.

Per family:
  1. Candidates (11_assign_families.py) covering >= 50% of their best reference hit.
  2. MAFFT --addfragments --keeplength onto the existing PLAZA representative
     alignment (rep.aln), so PLAZA columns, and the trimAl column mask computed on
     rep.aln, are unchanged; transcript fragments do not distort the alignment.
  3. Contamination filter: a new sequence is dropped if its most similar PLAZA
     sequence belongs to a different major lineage (fern / lycophyte /
     gymnosperm / monocot / eudicot) while PLAZA sequences of its own lineage
     are present in the family.
  4. Fragment filter: drop sequences with < 50% of trimmed columns filled.
  5. FastTree (LG+Gamma) optimizes branch lengths with the topology fixed
     (-nome -mllen -intree) on the expanded species tree (polytomies resolved once,
     seeded) pruned to the taxa present.
Output: data/work/gene_trees_x/<fam>/{aln.fa, topology.nwk, tree.nwk, qc.tsv, done}
Usage: python3 13_gene_trees_expanded.py [n|all] [workers]
"""
import pathlib, subprocess, sys, collections, random
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd, dendropy

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC = ROOT / "data/processed"
OLD = ROOT / "data/work/gene_trees"
WORK = ROOT / "data/work/gene_trees_x"; WORK.mkdir(parents=True, exist_ok=True)
CAND = ROOT / "data/work/candidates"
BIN = pathlib.Path("/home/user/tools/env/bin")
LIMIT = sys.argv[1] if len(sys.argv) > 1 else "all"
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 4
MIN_COV, MIN_FILL = 0.5, 0.5


def major(cls, order, sp=None):
    if cls == "Polypodiopsida": return "fern"
    if cls in ("Lycopodiopsida", "Isoetopsida"): return "lycophyte"
    if cls in ("Pinopsida", "Cycadopsida", "Ginkgoopsida", "Gnetopsida"): return "gymno"
    if cls == "Liliopsida": return "monocot"
    if cls == "Magnoliopsida":
        basal = {"Amborellales", "Nymphaeales", "Austrobaileyales", "Chloranthales", "Magnoliales", "Laurales",
                 "Piperales", "Canellales", "Ceratophyllales"}
        return "basal_angio" if order in basal else "eudicot"
    return "other"


def read_fasta(path):
    out, name = {}, None
    for line in open(path):
        line = line.rstrip("\n")
        if line.startswith(">"): name = line[1:].split()[0]; out[name] = []
        elif name: out[name].append(line)
    return {k: "".join(v) for k, v in out.items()}


def trim_mask(rep_aln, d):
    """Columns of rep.aln kept by trimAl -gappyout (via -colnumbering)."""
    r = subprocess.run([BIN / "trimal", "-in", rep_aln, "-out", d / "tmp.trim", "-gappyout", "-colnumbering"],
                       check=True, capture_output=True, text=True)
    (d / "tmp.trim").unlink(missing_ok=True)
    line = [l for l in r.stdout.splitlines() if l.startswith("#ColumnsMap")][0]
    return np.array([int(x) for x in line.split("\t", 1)[1].replace(",", " ").split()])


def run_family(args):
    fam, cands, groups, tree_path = args
    d = WORK / fam
    if (d / "done").exists(): return fam, "cached"
    d.mkdir(exist_ok=True)
    try:
        rep = OLD / fam / "rep.aln"
        ref = read_fasta(rep)
        with open(d / "new.fa", "w") as fh:
            for sp, seq in cands: fh.write(f">{sp}\n{seq}\n")
        with open(d / "add.aln", "w") as out:
            subprocess.run([BIN / "mafft", "--addfragments", d / "new.fa", "--keeplength", "--thread", "1", "--quiet", rep],
                           stdout=out, check=True)
        aln = read_fasta(d / "add.aln")
        cols = trim_mask(rep, d)
        names = list(aln); A = np.array([list(aln[n].upper()) for n in names])
        is_ref = np.array([n in ref for n in names])
        G = A != "-"
        qc = []
        keep = [i for i in range(len(names)) if is_ref[i]]
        ref_idx = np.where(is_ref)[0]
        ref_groups = {groups.get(names[i], "other") for i in ref_idx}
        for i in np.where(~is_ref)[0]:
            ov = G[i][None, :] & G[ref_idx]
            ident = ((A[i][None, :] == A[ref_idx]) & ov).sum(1) / np.maximum(ov.sum(1), 1)
            nn = names[ref_idx[ident.argmax()]]
            g_self, g_nn = groups.get(names[i], "other"), groups.get(nn, "other")
            fill = G[i, cols].mean()
            contam = g_self != g_nn and g_self in ref_groups and not ({g_self, g_nn} <= {"basal_angio", "monocot", "eudicot"} and "basal_angio" in (g_self, g_nn))
            status = "contam" if contam else ("fragment" if fill < MIN_FILL else "ok")
            qc.append((names[i], nn, round(float(ident.max()), 3), round(float(fill), 3), status))
            if status == "ok": keep.append(i)
        pd.DataFrame(qc, columns=["sp", "nearest_ref", "identity", "fill", "status"]).to_csv(d / "qc.tsv", sep="\t", index=False)
        T = A[np.ix_(keep, cols)]
        with open(d / "aln.fa", "w") as fh:
            for i, row in zip(keep, T): fh.write(f">{names[i]}\n{''.join(row)}\n")
        t = dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True, rooting="force-rooted")
        present = [names[i] for i in keep]
        t.retain_taxa_with_labels(present); t.suppress_unifurcations()
        (d / "topology.nwk").write_text(t.as_string(schema="newick").replace("[&R] ", "").strip())
        with open(d / "tree.nwk", "w") as out:
            subprocess.run([BIN / "FastTree", "-lg", "-gamma", "-nome", "-mllen", "-intree", d / "topology.nwk",
                            "-quiet", "-nopr", d / "aln.fa"], stdout=out, stderr=subprocess.DEVNULL, check=True)
        for f in ["new.fa", "add.aln"]: (d / f).unlink(missing_ok=True)
        (d / "done").write_text("ok")
        return fam, "ok"
    except Exception as e:
        return fam, f"error {e}"


if __name__ == "__main__":
    st = pd.read_csv(PROC / "species_table.tsv", sep="\t")
    # resolve Open Tree polytomies once (seeded) so all gene trees and the master tree share one topology
    tree_path = PROC / "species_tree_expanded_binary.nwk"
    if not tree_path.exists():
        t = dendropy.Tree.get(path=PROC / "species_tree_expanded.nwk", schema="newick", preserve_underscores=True,
                              rooting="force-rooted")
        t.resolve_polytomies(rng=random.Random(0))
        tree_path.write_text(t.as_string(schema="newick").replace("[&R] ", "").strip())
    tips = {l.taxon.label for l in dendropy.Tree.get(path=tree_path, schema="newick", preserve_underscores=True).leaf_node_iter()}
    groups = {sp: major(c, o) for sp, c, o in zip(st.sp, st["class"], st.order)}
    fams = (PROC / "rer_families.txt").read_text().split()
    if LIMIT != "all": fams = fams[:int(LIMIT)]
    todo = [f for f in fams if not (WORK / f / "done").exists()]
    print(f"{len(fams)} families, {len(todo)} to run", flush=True)
    cands = collections.defaultdict(list)
    for p in CAND.glob("*.tsv"):
        sp = p.stem
        if sp not in tips: continue
        c = pd.read_csv(p, sep="\t")
        c = c[(c.alen / c.slen >= MIN_COV) & c.family.isin(todo)]
        for f, s in zip(c.family, c.seq): cands[f].append((sp, s))
    jobs = [(f, cands.get(f, []), groups, tree_path) for f in todo]
    status = collections.Counter()
    with ProcessPoolExecutor(NW) as ex:
        for i, (fam, s) in enumerate(ex.map(run_family, jobs, chunksize=1)):
            status[s.split()[0]] += 1
            if s.startswith("error"): print(fam, s[:300], flush=True)
            if i % 50 == 0: print(i, dict(status), flush=True)
    print("finished", dict(status))
