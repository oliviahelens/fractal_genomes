"""Build a unified species x gene-family copy-number matrix.

Backbone: PLAZA dicots 5.0 homologous families (HOM05D...).
Other instances (monocots 5.0, basal 1.0 [ferns], gymno 1.0) define their own
families. Each foreign family is mapped to the backbone family that holds the
most of its genes from species shared between the two instances (identical
gene IDs). Foreign families with no anchor genes are kept as new families.
"""
import gzip, collections, pathlib
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
RAW = ROOT / "data/raw/plaza"
OUT = ROOT / "data/processed"; OUT.mkdir(parents=True, exist_ok=True)


def read_fam(path):
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#") or line.startswith('"gf_id"'):
                continue
            sep = "\t" if "\t" in line else ";"
            gf, sp, gene = [x.strip().strip('"') for x in line.rstrip("\n").split(sep)[:3]]
            rows.append((gf, sp, gene))
    return pd.DataFrame(rows, columns=["gf", "sp", "gene"])


backbone = read_fam(RAW / "dicots_HOMFAM.csv.gz")
gene2bb = dict(zip(backbone.gene, backbone.gf))
bb_species = set(backbone.sp)

# which species to take from which instance (first instance that has it wins)
instances = [
    ("monocots", RAW / "plaza_public_monocots_05_HOMFAM.csv.gz"),
    ("basal", RAW / "plaza_basal_01_HOMFAM.csv.gz"),
    ("gymno", RAW / "plaza_gymno_01_HOMFAM.csv.gz"),
]
assignments = [backbone.assign(instance="dicots")]
taken = set(bb_species)
log = []
for name, path in instances:
    df = read_fam(path)
    anchors = df[df.gene.isin(gene2bb)]
    # foreign family -> backbone family by majority vote of anchor genes
    vote = (anchors.assign(bb=anchors.gene.map(gene2bb))
            .groupby(["gf", "bb"]).size().reset_index(name="n")
            .sort_values("n", ascending=False).drop_duplicates("gf"))
    fmap = dict(zip(vote.gf, vote.bb))
    new_sp = sorted(set(df.sp) - taken)
    sub = df[df.sp.isin(new_sp)].copy()
    sub["gf"] = sub.gf.map(lambda g: fmap.get(g, f"{name}:{g}"))
    sub["instance"] = name
    assignments.append(sub)
    taken |= set(new_sp)
    mapped = sub.gf.str.startswith("HOM05D").mean() if len(sub) else 0
    log.append(f"{name}: added {len(new_sp)} species {new_sp}; "
               f"{len(fmap)}/{df.gf.nunique()} families mapped; "
               f"{mapped:.1%} of new-species genes in backbone families")

allg = pd.concat(assignments, ignore_index=True)
allg[["gene", "sp", "gf"]].to_csv(OUT / "gene_family_assignments.tsv.gz", sep="\t", index=False, compression="gzip")
counts = allg.groupby(["gf", "sp"]).size().rename("n").reset_index()
counts.to_csv(OUT / "family_counts_long.tsv.gz", sep="\t", index=False, compression="gzip")
# analysis matrix: backbone families present in >= 3 species (keeps memory small)
nsp = counts.groupby("gf").sp.nunique()
keep = nsp[(nsp >= 3) & nsp.index.str.startswith("HOM05D")].index
mat = counts[counts.gf.isin(keep)].pivot(index="sp", columns="gf", values="n").fillna(0).astype("int32")
mat.insert(0, "_total_genes", allg.groupby("sp").size().reindex(mat.index))
mat.to_csv(OUT / "family_copy_number.tsv.gz", sep="\t", compression="gzip")
allg[["sp", "instance"]].drop_duplicates().to_csv(OUT / "species_instance.tsv", sep="\t", index=False)
(OUT / "family_matrix_log.txt").write_text("\n".join(log) + f"\nmatrix: {mat.shape}\n")
print("\n".join(log)); print("matrix", mat.shape)
