"""Leaf/branching regulator families (approach A): collect every copy, define subclades.

Stage `collect`:
  * reference = all members of the 28 regulator families (data/regulators/anchors.tsv)
    in the 30 reference PLAZA genomes; full reference = their complete proteomes;
  * every species of the expanded set: proteins hitting a regulator reference (DIAMOND,
    e < 1e-10, >= 50% of the shorter sequence aligned) are kept only if their best hit in
    the complete reference proteomes belongs to the same regulator family (reciprocal check:
    removes proteins that merely share a domain); PLAZA species use their PLAZA family
    assignments directly;
  -> data/work/regulators/fwd_all.tsv.gz (sp, prot, family, best_rid, bits, seq)
Stage `subclades`:
  * per family, a FastTree tree of the reference members; each reference gene is labelled
    with the subclade of its nearest Arabidopsis gene (patristic) if that gene is an anchor;
  * each member inherits the label of its best reference hit;
  -> data/work/regulators/ref_labels.tsv, ref_trees/<family>.nwk
Stage `verify`: members whose best hit is an anchor-subclade reference gene, reciprocal check
  against the complete reference proteomes -> data/work/regulators/members_subclade.tsv
Usage: python3 17_regulators.py collect|subclades|verify [workers]
"""
import gzip, pathlib, subprocess, sys, collections
from concurrent.futures import ProcessPoolExecutor
import numpy as np, pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC = ROOT / "data/processed"
WORK = ROOT / "data/work/regulators"; WORK.mkdir(parents=True, exist_ok=True)
BIN = pathlib.Path("/home/user/tools/env/bin")
STAGE = sys.argv[1] if len(sys.argv) > 1 else "collect"
NW = int(sys.argv[2]) if len(sys.argv) > 2 else 4
REF_SP = ["ath", "ptr", "gma", "vvi", "egr", "csa", "sly", "oeu", "han", "dca", "bvu", "nnu", "aox", "pso",
          "atr", "cde", "mbi", "osa", "zma", "mac", "egu", "aco", "zosmarina", "gbi", "pab", "cmi", "gmo",
          "smo", "cri", "afi"]
anchors = pd.read_csv(ROOT / "data/regulators/anchors.tsv", sep="\t")
ga = pd.read_csv(PROC / "gene_family_assignments.tsv.gz", sep="\t")
ath = ga[ga.sp == "ath"].assign(g=lambda d: d.gene.str.upper().str.replace(r"\..*", "", regex=True))
anchors["family"] = anchors.ath_gene.map(dict(zip(ath.g, ath.gf)))
FAMS = sorted(anchors.family.unique())


def read_fasta_gz(path, key):
    out, name, buf = {}, None, []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name: out[name] = "".join(buf).rstrip("*")
                name, buf = key(line[1:].strip()), []
            else:
                buf.append(line.strip())
    if name: out[name] = "".join(buf).rstrip("*")
    return out


plaza_key = lambda h: h.split("|")[1].strip() if "|" in h else h.split()[0]


def plaza_proteome(sp):
    return read_fasta_gz(next((ROOT / "data/raw/proteomes").glob(f"{sp}.*gz")), plaza_key)


def build_refs():
    if (WORK / "full_ref.dmnd").exists():
        return pd.read_csv(WORK / "reg_ref.tsv", sep="\t"), pd.read_csv(WORK / "full_ref.tsv", sep="\t")
    fam_of = {(s, g): f for s, g, f in zip(ga.sp, ga.gene, ga.gf)}
    reg_rows, full_rows = [], []
    with open(WORK / "reg_ref.fa", "w") as fr, open(WORK / "full_ref.fa", "w") as ff:
        for sp in REF_SP:
            seqs = plaza_proteome(sp)
            for g, s in seqs.items():
                f = fam_of.get((sp, g), "")
                fid = f"f{len(full_rows)}"; full_rows.append((fid, sp, g, f)); ff.write(f">{fid}\n{s}\n")
                if f in FAMS:
                    rid = f"r{len(reg_rows)}"; reg_rows.append((rid, sp, f, g)); fr.write(f">{rid}\n{s}\n")
    reg = pd.DataFrame(reg_rows, columns=["rid", "sp", "family", "gene"]); reg.to_csv(WORK / "reg_ref.tsv", sep="\t", index=False)
    full = pd.DataFrame(full_rows, columns=["fid", "sp", "gene", "family"]); full.to_csv(WORK / "full_ref.tsv", sep="\t", index=False)
    for name in ["reg_ref", "full_ref"]:
        subprocess.run([BIN / "diamond", "makedb", "--in", WORK / f"{name}.fa", "-d", WORK / name, "--quiet"], check=True)
    return reg, full


def forward(args):
    """Proteins of one species with a regulator hit: best reference hit per protein."""
    sp, path, kind = args
    out = WORK / "fwd" / f"{sp}.tsv"
    if out.exists(): return sp
    seqs = read_fasta_gz(path, plaza_key if kind == "plaza" else (lambda h: h.split()[0]))
    q = WORK / "fwd" / f"{sp}.fa"
    with open(q, "w") as fh:
        for k, s in seqs.items():
            if len(s) >= 30: fh.write(f">{k}\n{s}\n")
    m8 = WORK / "fwd" / f"{sp}.m8"
    subprocess.run([BIN / "diamond", "blastp", "-q", q, "-d", WORK / "reg_ref", "-o", m8, "--evalue", "1e-10",
                    "--max-target-seqs", "1", "--threads", "1", "--quiet",
                    "--outfmt", "6", "qseqid", "sseqid", "bitscore", "qlen", "slen", "length"], check=True)
    h = pd.read_csv(m8, sep="\t", header=None, names=["prot", "best_rid", "bits", "qlen", "slen", "alen"])
    h = h.sort_values("bits", ascending=False).drop_duplicates("prot")
    h = h[h.alen >= 0.5 * np.minimum(h.qlen, h.slen)]
    h["seq"] = h.prot.map(seqs); h.insert(0, "sp", sp)
    h.to_csv(out, sep="\t", index=False)
    q.unlink(); m8.unlink()
    return sp


if STAGE == "collect":
    reg, full = build_refs()
    rid_fam = dict(zip(reg.rid, reg.family))
    (WORK / "fwd").mkdir(exist_ok=True)
    sp_all = pd.read_csv(ROOT / "results/rer_expanded/rer_cv_predictions.tsv", sep="\t").sp
    st = pd.read_csv(PROC / "species_table.tsv", sep="\t").set_index("sp").loc[sp_all]
    jobs = []
    for sp, row in st.iterrows():
        if row.source == "plaza":
            jobs.append((sp, next((ROOT / "data/raw/proteomes").glob(f"{sp}.*gz")), "plaza"))
        elif row.source == "1kp":
            jobs.append((sp, ROOT / f"data/raw/onekp/{sp}.faa.gz", "1kp"))
        else:
            jobs.append((sp, ROOT / f"data/raw/ncbi/{row.accession}.faa.gz", "ncbi"))
    with ProcessPoolExecutor(NW) as ex:
        for i, sp in enumerate(ex.map(forward, jobs, chunksize=4)):
            if i % 100 == 0: print("forward", i, sp, flush=True)
    fwd = pd.concat([pd.read_csv(WORK / "fwd" / f"{sp}.tsv", sep="\t") for sp, _, _ in jobs], ignore_index=True)
    fwd["family"] = fwd.best_rid.map(rid_fam)
    fwd.to_csv(WORK / "fwd_all.tsv.gz", sep="\t", index=False)
    print("forward hits", len(fwd))

elif STAGE == "subclades":
    import dendropy
    reg = pd.read_csv(WORK / "reg_ref.tsv", sep="\t")
    (WORK / "ref_trees").mkdir(exist_ok=True)
    seqs = {}
    cur = None
    for line in open(WORK / "reg_ref.fa"):
        if line.startswith(">"): cur = line[1:].strip()
        else: seqs[cur] = line.strip()
    anchor_sub = dict(zip(anchors.ath_gene, anchors.subclade))
    label = {}
    for fam in FAMS:
        sub = reg[reg.family == fam]
        tre = WORK / "ref_trees" / f"{fam}.nwk"
        if not tre.exists():
            fa = WORK / "ref_trees" / f"{fam}.fa"
            with open(fa, "w") as fh:
                for r in sub.rid: fh.write(f">{r}\n{seqs[r]}\n")
            with open(WORK / "ref_trees" / f"{fam}.aln", "w") as out:
                subprocess.run([BIN / "mafft", "--auto", "--thread", str(NW), "--quiet", fa], stdout=out, check=True)
            with open(tre, "w") as out:
                subprocess.run([BIN / "FastTreeMP", "-lg", "-gamma", "-quiet", "-nopr", WORK / "ref_trees" / f"{fam}.aln"],
                               stdout=out, stderr=subprocess.DEVNULL, check=True)
        t = dendropy.Tree.get(path=tre, schema="newick")
        pdm = t.phylogenetic_distance_matrix()
        leaves = {l.taxon.label: l.taxon for l in t.leaf_node_iter()}
        ath_r = [r for r, s, g in zip(sub.rid, sub.sp, sub.gene) if s == "ath"]
        ath_gene = {r: g.upper().split(".")[0] for r, s, g in zip(sub.rid, sub.sp, sub.gene) if s == "ath"}
        for r in sub.rid:
            if r not in leaves: continue
            d = [(pdm.patristic_distance(leaves[r], leaves[a]), a) for a in ath_r if a in leaves]
            nearest = min(d)[1]
            label[r] = anchor_sub.get(ath_gene[nearest], "other")
        print(fam, collections.Counter(label[r] for r in sub.rid if r in label).most_common(6), flush=True)
    pd.DataFrame({"rid": list(label), "subclade": list(label.values())}).to_csv(WORK / "ref_labels.tsv", sep="\t", index=False)

elif STAGE == "verify":
    # keep proteins whose best regulator hit is a reference gene of an anchor subclade, then require
    # their best hit among the complete reference proteomes to be in the same family (chunked search)
    reg, full = build_refs()
    lab = pd.read_csv(WORK / "ref_labels.tsv", sep="\t")
    lab = dict(zip(lab.rid, lab.subclade))
    fwd = pd.read_csv(WORK / "fwd_all.tsv.gz", sep="\t")
    fwd["subclade"] = fwd.best_rid.map(lab).fillna("other")
    fwd = fwd[fwd.subclade != "other"].reset_index(drop=True)
    st = pd.read_csv(PROC / "species_table.tsv", sep="\t").set_index("sp")
    is_pl = fwd.sp.map(st.source) == "plaza"
    pl_fam = {(s_, g): f for s_, g, f in zip(ga.sp, ga.gene, ga.gf)}
    keep = np.zeros(len(fwd), bool)
    keep[is_pl.values] = [pl_fam.get((s_, p_)) == f for s_, p_, f in zip(fwd.sp[is_pl], fwd.prot[is_pl], fwd.family[is_pl])]
    fid_fam = dict(zip(full.fid, full.family))
    oth = np.flatnonzero(~is_pl.values)
    print("candidates", len(fwd), "plaza", is_pl.sum(), "to verify", len(oth), flush=True)
    (WORK / "tmp").mkdir(exist_ok=True)
    for c0 in range(0, len(oth), 50000):
        idx = oth[c0:c0 + 50000]; qfa = WORK / "tmp" / "rq.fa"; m8 = WORK / "tmp" / f"rv{c0}.m8"
        if not m8.exists():
            with open(qfa, "w") as fh:
                for i in idx: fh.write(f">q{i}\n{fwd.seq[i]}\n")
            subprocess.run([BIN / "diamond", "blastp", "-q", qfa, "-d", WORK / "full_ref", "-o", m8.with_suffix(".part"),
                            "--evalue", "1e-10", "--max-target-seqs", "1", "--max-hsps", "1", "--threads", str(NW),
                            "--quiet", "-t", WORK / "tmp", "--outfmt", "6", "qseqid", "sseqid", "bitscore"], check=True)
            m8.with_suffix(".part").rename(m8)
        rv = pd.read_csv(m8, sep="\t", header=None, names=["q", "fid", "bits"]).sort_values("bits", ascending=False).drop_duplicates("q")
        best = dict(zip(rv.q, rv.fid.map(fid_fam)))
        for i in idx: keep[i] = best.get(f"q{i}") == fwd.family[i]
        print("verified", c0 + len(idx), flush=True)
    members = fwd[keep]
    members.to_csv(WORK / "members_subclade.tsv", sep="\t", index=False)
    print("kept", len(members), "of", len(fwd))
    print(members.groupby("subclade").size().to_string())
