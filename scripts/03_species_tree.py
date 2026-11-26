"""Fetch a species tree for all species from the Open Tree of Life synthetic tree.

Open Tree topology has no branch lengths; Grafen (1989) branch lengths are
applied downstream. Output: data/processed/species_tree.nwk with tips = PLAZA codes.
"""
import json, pathlib, re
import pandas as pd, requests

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC = ROOT / "data/processed"
API = "https://api.opentreeoflife.org/v3"

names = pd.read_csv(ROOT / "data/raw/plaza/species_names.tsv", sep="\t", header=None,
                    names=["sp", "name", "taxid"])
extra = pd.DataFrame([  # species from instances without SpeciesInformation files
    ("afi", "Azolla filiculoides"), ("cri", "Ceratopteris richardii"), ("scu", "Salvinia cucullata"),
    ("cmi", "Cycas micholitzii"), ("gbi", "Ginkgo biloba"), ("gmo", "Gnetum montanum"),
    ("pab", "Picea abies"), ("pgl", "Picea glauca"), ("pme", "Pseudotsuga menziesii"),
    ("ppi", "Pinus pinaster"), ("psi", "Picea sitchensis"), ("psy", "Pinus sylvestris"),
    ("pta", "Pinus taeda"), ("tba", "Taxus baccata")], columns=["sp", "name"])
sp = pd.concat([names[["sp", "name"]], extra]).drop_duplicates("sp")
sp["query"] = sp.name.str.replace(r"\s+(L\.|ssp\..*|var\..*|Kitaake)$", "", regex=True).str.strip()
sp = sp[~sp.sp.isin(["osaindica", "osak", "zma-b104"])]  # duplicate genomes of one species

r = requests.post(f"{API}/tnrs/match_names", json={"names": sp["query"].tolist(), "do_approximate_matching": False}, timeout=120)
r.raise_for_status()
ott = {}
for res in r.json()["results"]:
    if res["matches"]:
        ott[res["name"]] = res["matches"][0]["taxon"]["ott_id"]
sp["ott_id"] = sp["query"].map(ott)
missing = sp[sp.ott_id.isna()]
print("unmatched:", missing[["sp", "query"]].values.tolist())
sp = sp.dropna(subset=["ott_id"]).astype({"ott_id": int})

r = requests.post(f"{API}/tree_of_life/induced_subtree", json={"ott_ids": sp.ott_id.tolist(), "label_format": "id"}, timeout=120)
if r.status_code != 200:  # some ott ids may be pruned from synthesis; drop and retry
    body = r.json(); print(body.get("message", "")[:500])
    bad = {int(x) for x in re.findall(r"ott(\d+)", json.dumps(body.get("unknown", {})))}
    print("dropped from synthesis:", sp[sp.ott_id.isin(bad)].sp.tolist())
    sp_all = sp.copy(); sp = sp[~sp.ott_id.isin(bad)]
    r = requests.post(f"{API}/tree_of_life/induced_subtree", json={"ott_ids": sp.ott_id.tolist(), "label_format": "id"}, timeout=120)
r.raise_for_status()
res = r.json(); nwk = res["newick"]
lookup = dict(zip("ott" + sp.ott_id.astype(str), sp.sp))
for ottid, mrca in res.get("broken", {}).items():  # taxa not monophyletic in synthesis
    lookup[mrca] = lookup[ottid]
nwk = re.sub(r"([(,])(mrcaott\d+ott\d+|ott\d+)", lambda m: m.group(1) + lookup[m.group(2)], nwk)
nwk = re.sub(r"\)[^,;()]+", ")", nwk)  # drop internal node labels
import dendropy
t = dendropy.Tree.get(data=nwk, schema="newick", preserve_underscores=True)
t.suppress_unifurcations()
nwk = t.as_string(schema="newick").strip().replace("[&R] ", "")
# Salvia bowleyana (sbo) is absent from synthesis; graft it into Lamiales next to
# Striga + Erythranthe (Lamiaceae/Orobanchaceae/Phrymaceae clade).
if "sbo" not in nwk and "(sas,egut)" in nwk:
    nwk = nwk.replace("(sas,egut)", "((sas,egut),sbo)")
(PROC / "species_tree.nwk").write_text(nwk)
sp = pd.concat([sp, sp_all[sp_all.sp == "sbo"]]) if "sbo" in nwk else sp
sp.to_csv(PROC / "species_names.tsv", sep="\t", index=False)
print(nwk[:600]); print(len(sp), "species")
