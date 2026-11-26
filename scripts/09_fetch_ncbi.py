"""Download NCBI annotated reference genomes (vascular plants) for species not in PLAZA.

One assembly per species (NCBI 'reference' flag; RefSeq preferred over GenBank).
Writes data/raw/ncbi/<accession>.faa.gz (annotated proteins) and genomes.tsv
(accession, species, taxid, annotation provider, protein-coding genes).
Usage: python3 09_fetch_ncbi.py [n_workers]
"""
import io, sys, zipfile, gzip, pathlib, time, requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/ncbi"; OUT.mkdir(parents=True, exist_ok=True)
API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
NW = int(sys.argv[1]) if len(sys.argv) > 1 else 3


def get(url, **kw):
    for i in range(6):
        try:
            r = requests.get(url, timeout=600, **kw); r.raise_for_status(); return r
        except Exception:
            if i == 5: raise
            time.sleep(10 * (i + 1))


def reports():
    tok, out = None, []
    while True:
        p = {"filters.has_annotation": "true", "filters.reference_only": "true", "page_size": 500}
        if tok: p["page_token"] = tok
        d = get(f"{API}/genome/taxon/Tracheophyta/dataset_report", params=p).json()
        out += d.get("reports", []); tok = d.get("next_page_token")
        if not tok: return out


def fetch(acc):
    dest = OUT / f"{acc}.faa.gz"
    if dest.exists(): return acc, "cached"
    try:
        z = zipfile.ZipFile(io.BytesIO(get(f"{API}/genome/accession/{acc}/download",
                                           params={"include_annotation_type": "PROT_FASTA"}).content))
        name = next((n for n in z.namelist() if n.endswith("protein.faa")), None)
        if name is None: return acc, "no proteins"
        tmp = dest.with_suffix(".tmp")
        with gzip.open(tmp, "wb") as fh: fh.write(z.read(name))
        tmp.rename(dest); return acc, "ok"
    except Exception as e:
        return acc, f"error {e}"


if __name__ == "__main__":
    rows = []
    for r in reports():
        a = r.get("annotation_info", {})
        name = r["organism"]["organism_name"]
        rows.append({"accession": r["accession"], "organism": name, "species": " ".join(name.split()[:2]),
                     "taxid": r["organism"]["tax_id"], "provider": a.get("provider", ""),
                     "genes": a.get("stats", {}).get("gene_counts", {}).get("protein_coding", 0)})
    g = pd.DataFrame(rows)
    plaza = set(pd.read_csv(ROOT / "data/raw/plaza/species_names.tsv", sep="\t", header=None)[1])
    g["in_plaza"] = g.species.isin(plaza)
    g["refseq"] = g.accession.str.startswith("GCF_")
    g = g[~g.in_plaza & (g.genes > 5000)].sort_values(["species", "refseq", "genes"], ascending=[True, False, False])
    g = g.drop_duplicates("species")
    g.to_csv(OUT / "genomes.tsv", sep="\t", index=False)
    print(f"{len(g)} new species ({g.refseq.sum()} RefSeq)", flush=True)
    st = {}
    with ThreadPoolExecutor(NW) as ex:
        for i, (acc, s) in enumerate(ex.map(fetch, g.accession)):
            st[s.split()[0]] = st.get(s.split()[0], 0) + 1
            if not s.startswith(("ok", "cached")): print(acc, s[:200], flush=True)
            if i % 50 == 0: print(i, st, flush=True)
    print("finished", st)
