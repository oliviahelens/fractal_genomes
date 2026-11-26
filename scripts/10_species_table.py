"""Merge PLAZA genomes, NCBI genomes and 1KP transcriptomes into one species table.

One entry per species. Source priority: PLAZA genome > NCBI RefSeq > NCBI GenBank
annotation > 1KP transcriptome (best BUSCO completeness among that species'
samples). 1KP samples are dropped if the 1KP consortium flagged them as
multi-species, misplaced in the tree, or contaminated, and if they are not
vascular plants. Adds NCBI taxonomy (family, order) for every species.
Output: data/processed/species_table.tsv
"""
import pathlib, re, time, requests
import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
META = ROOT / "data/raw/onekp/meta"
VASC = {"Eudicots", "Monocots", "Gymnos", "Monilophytes", "Magnoliids", "Lycophytes", "ANAGrade",
        "Chloranthales", "Ceratophyllales"}


def binomial(s):
    s = re.sub(r"\b(cf|aff)\.\s*", "", str(s)).replace("_", " ").replace("×", "x").strip()
    w = s.split()
    return " ".join(w[:3] if len(w) > 2 and w[1] == "x" else w[:2])


def flagged():
    bad, sec = set(), None
    keep_secs = ("No Family-level", "No SSU", "Unconfirmed Source (No Worrisome", "Status Changed")
    for line in open(META / "contamination.txt"):
        line = line.strip()
        if line.startswith("#"):
            if not line.lstrip("# ").isupper() or len(line.lstrip("# ")) != 4:
                sec = line.lstrip("# "); continue
        code = line.lstrip("# ").strip()
        if len(code) == 4 and code.isupper() and not (sec or "").startswith(keep_secs):
            bad.add(code)
    return bad


def taxonomy(names):
    out = {}
    for i in range(0, len(names), 200):
        batch = names[i:i + 200]
        for k in range(5):
            try:
                r = requests.post("https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/dataset_report",
                                  json={"taxons": batch}, timeout=120).json(); break
            except Exception:
                time.sleep(5 * (k + 1))
        for rep in r.get("reports", []):
            t = rep.get("taxonomy", {}); c = t.get("classification", {})
            for q in rep.get("query", []):
                out[q] = {"family": c.get("family", {}).get("name", ""), "order": c.get("order", {}).get("name", ""),
                          "class": c.get("class", {}).get("name", ""), "taxid": t.get("tax_id", "")}
    return out


plaza = pd.read_csv(ROOT / "data/processed/species_names.tsv", sep="\t").rename(columns={"name": "species"})
plaza = plaza[["sp", "species"]].assign(taxid="")
plaza = plaza.assign(source="plaza", data_type="genome", completeness=float("nan"))

ncbi = pd.read_csv(ROOT / "data/raw/ncbi/genomes.tsv", sep="\t")
ncbi = pd.DataFrame({"sp": "N" + ncbi.accession.str.replace(r"\..*", "", regex=True).str.replace("_", ""),
                     "species": ncbi.organism.map(binomial), "taxid": ncbi.taxid,
                     "source": ["ncbi_refseq" if r else "ncbi_genbank" for r in ncbi.refseq],
                     "data_type": "genome", "completeness": float("nan"), "accession": ncbi.accession})

ann = pd.read_csv(META / "annotations.csv")
busco = pd.read_csv(META / "CEGMA_BUSCO_100-missing.tsv", sep="\t").rename(columns={"V1": "code", "V2": "busco"})
samples = pd.read_csv(ROOT / "data/raw/onekp/samples.tsv", sep="\t")
ok = ann[ann["Very Brief Classification"].isin(VASC) & ~ann.Code.isin(flagged()) & ann.Code.isin(samples.code)]
onekp = pd.DataFrame({"sp": ok.Code.values, "species": ok.Species.map(binomial).values, "taxid": "",
                      "source": "1kp", "data_type": "transcriptome",
                      "onekp_group": ok["Very Brief Classification"].values})
onekp["completeness"] = onekp.sp.map(dict(zip(busco.code, busco.busco)))
onekp = onekp[~onekp.species.str.contains(r"\bsp\.?$|\bsp\.? |unknown", case=False, regex=True)]
onekp = onekp.sort_values("completeness", ascending=False).drop_duplicates("species")

allsp = pd.concat([plaza, ncbi, onekp], ignore_index=True)
allsp["species"] = allsp.species.map(binomial)
prio = {"plaza": 0, "ncbi_refseq": 1, "ncbi_genbank": 2, "1kp": 3}
allsp["prio"] = allsp.source.map(prio)
allsp = allsp.sort_values(["prio"]).drop_duplicates("species").drop(columns="prio")
tx = taxonomy(sorted(allsp.species.unique()))
for col in ["family", "order", "class"]:
    allsp[col] = allsp.species.map(lambda s: tx.get(s, {}).get(col, ""))
allsp["taxid"] = [t if str(t) not in ("", "nan") else tx.get(s, {}).get("taxid", "") for t, s in zip(allsp.taxid, allsp.species)]
allsp.to_csv(ROOT / "data/processed/species_table.tsv", sep="\t", index=False)
print(allsp.source.value_counts().to_string())
print("no NCBI taxonomy:", (allsp.family == "").sum())
