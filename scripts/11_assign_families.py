"""Assign NCBI and 1KP proteins to the 3,989 RER gene families.

Reference: all members of those families from 30 PLAZA genomes chosen to span
the tree (ferns, lycophyte, gymnosperms, ANA, monocots, eudicots). Each new
species' proteins are searched with DIAMOND (e < 1e-10); a protein is assigned
to the family of its best hit. Per species and family the protein with the
highest bitscore is kept as the candidate (isoforms/paralogs/fragments: the
best-matching one wins). Completeness = fraction of families with a candidate.
Outputs: data/work/candidates/<sp>.tsv (family, protein id, bitscore, pident,
qlen, slen, ref species, seq) and data/processed/completeness.tsv.
Usage: python3 11_assign_families.py [threads] [species...]
"""
import gzip, pathlib, subprocess, sys, collections
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC = ROOT / "data/processed"
OUT = ROOT / "data/work/candidates"; OUT.mkdir(parents=True, exist_ok=True)
DB = ROOT / "data/work/ref_families"
BIN = pathlib.Path("/home/user/tools/env/bin")
TH = sys.argv[1] if len(sys.argv) > 1 else "4"
REF_SP = ["ath", "ptr", "gma", "vvi", "egr", "csa", "sly", "oeu", "han", "dca", "bvu", "nnu", "aox", "pso",
          "atr", "cde", "mbi", "osa", "zma", "mac", "egu", "aco", "zosmarina", "gbi", "pab", "cmi", "gmo",
          "smo", "cri", "afi"]


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


def build_ref():
    if (DB.with_suffix(".dmnd")).exists():
        return pd.read_csv(DB.with_suffix(".tsv"), sep="\t")
    fams = set((PROC / "rer_families.txt").read_text().split())
    ga = pd.read_csv(PROC / "gene_family_assignments.tsv.gz", sep="\t")
    ga = ga[ga.gf.isin(fams) & ga.sp.isin(REF_SP)]
    rows = []
    with open(DB.with_suffix(".fa"), "w") as fh:
        for sp, sub in ga.groupby("sp"):
            f = next((ROOT / "data/raw/proteomes").glob(f"{sp}.*gz"))
            seqs = read_fasta_gz(f, lambda h: h.split("|")[1].strip() if "|" in h else h.split()[0])
            for g, gf in zip(sub.gene, sub.gf):
                if g in seqs:
                    rid = f"r{len(rows)}"; rows.append((rid, sp, gf, g))
                    fh.write(f">{rid}\n{seqs[g]}\n")
    ref = pd.DataFrame(rows, columns=["rid", "sp", "family", "gene"])
    ref.to_csv(DB.with_suffix(".tsv"), sep="\t", index=False)
    subprocess.run([BIN / "diamond", "makedb", "--in", DB.with_suffix(".fa"), "-d", DB, "--quiet"], check=True)
    return ref


def proteome_path(row):
    if row.source == "1kp":
        return ROOT / f"data/raw/onekp/{row.sp}.faa.gz"
    return ROOT / f"data/raw/ncbi/{row.accession}.faa.gz"


def run_species(row, ref_fam, ref_sp):
    out = OUT / f"{row.sp}.tsv"
    if out.exists():
        return
    src = proteome_path(row)
    if not src.exists():
        return
    seqs = read_fasta_gz(src, lambda h: h.split()[0])
    q = OUT / f"{row.sp}.query.fa"
    if not any(len(s) >= 30 for s in seqs.values()):  # empty 1KP archive
        pd.DataFrame(columns=["family", "q", "bits", "pident", "qlen", "slen", "alen", "ref_sp", "seq"]).to_csv(out, sep="\t", index=False)
        return
    with open(q, "w") as fh:
        for k, s in seqs.items():
            if len(s) >= 30: fh.write(f">{k}\n{s}\n")
    hits = OUT / f"{row.sp}.m8"
    subprocess.run([BIN / "diamond", "blastp", "-q", q, "-d", DB, "-o", hits, "--evalue", "1e-10",
                    "--max-target-seqs", "1", "--threads", TH, "--quiet",
                    "--outfmt", "6", "qseqid", "sseqid", "pident", "bitscore", "qlen", "slen", "length"], check=True)
    h = pd.read_csv(hits, sep="\t", header=None, names=["q", "rid", "pident", "bits", "qlen", "slen", "alen"])
    h = h.sort_values("bits", ascending=False).drop_duplicates("q")
    h["family"] = h.rid.map(ref_fam); h["ref_sp"] = h.rid.map(ref_sp)
    h = h.drop_duplicates("family")  # best protein per family
    h["seq"] = h.q.map(seqs)
    h[["family", "q", "bits", "pident", "qlen", "slen", "alen", "ref_sp", "seq"]].to_csv(out, sep="\t", index=False)
    q.unlink(); hits.unlink()


if __name__ == "__main__":
    ref = build_ref()
    ref_fam, ref_sp = dict(zip(ref.rid, ref.family)), dict(zip(ref.rid, ref.sp))
    st = pd.read_csv(PROC / "species_table.tsv", sep="\t")
    st = st[st.source != "plaza"]
    if len(sys.argv) > 2:
        st = st[st.sp.isin(sys.argv[2:])]
    nfam = len((PROC / "rer_families.txt").read_text().split())
    for i, row in enumerate(st.itertuples()):
        run_species(row, ref_fam, ref_sp)
        if i % 25 == 0: print(i, row.sp, flush=True)
    comp = [(p.stem, sum(1 for _ in open(p)) - 1) for p in OUT.glob("*.tsv")]
    pd.DataFrame([(s, n, n / nfam) for s, n in comp], columns=["sp", "families", "completeness"]).to_csv(
        PROC / "completeness.tsv", sep="\t", index=False)
    print("finished", len(comp))
