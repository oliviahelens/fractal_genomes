"""Open Tree of Life topology for the expanded species set (PLAZA + NCBI + 1KP).

Species unmatched by TNRS or absent from the synthetic tree are grafted as
sister to their congeners (MRCA of all congeners) when any congener is present;
otherwise dropped. Output: data/processed/species_tree_expanded.nwk (tips = sp
codes) and species_tree_expanded_log.tsv.
"""
import json, pathlib, re
import pandas as pd, requests, dendropy

ROOT = pathlib.Path(__file__).resolve().parents[1]
PROC = ROOT / "data/processed"
API = "https://api.opentreeoflife.org/v3"

st = pd.read_csv(PROC / "species_table.tsv", sep="\t")
ott = {}
names = st.species.unique().tolist()
for i in range(0, len(names), 500):
    r = requests.post(f"{API}/tnrs/match_names", json={"names": names[i:i + 500], "do_approximate_matching": False,
                                                        "context_name": "Land plants"}, timeout=300)
    r.raise_for_status()
    for res in r.json()["results"]:
        if res["matches"]:
            ott[res["name"]] = res["matches"][0]["taxon"]["ott_id"]
st["ott_id"] = st.species.map(ott)
st["status"] = ["tnrs_unmatched" if pd.isna(o) else "" for o in st.ott_id]
cand = st.dropna(subset=["ott_id"]).astype({"ott_id": int}).drop_duplicates("ott_id")
st.loc[st.ott_id.notna() & ~st.sp.isin(cand.sp), "status"] = "duplicate_ott"

ids = cand.ott_id.tolist()
r = requests.post(f"{API}/tree_of_life/induced_subtree", json={"ott_ids": ids, "label_format": "id"}, timeout=600)
if r.status_code != 200:
    bad = {int(x) for x in re.findall(r"ott(\d+)", json.dumps(r.json().get("unknown", {})))}
    st.loc[st.ott_id.isin(bad), "status"] = "not_in_synthesis"
    cand = cand[~cand.ott_id.isin(bad)]
    r = requests.post(f"{API}/tree_of_life/induced_subtree", json={"ott_ids": cand.ott_id.tolist(), "label_format": "id"}, timeout=600)
r.raise_for_status()
res = r.json(); nwk = res["newick"]
lookup = dict(zip("ott" + cand.ott_id.astype(str), cand.sp))
for ottid, mrca in res.get("broken", {}).items():
    if ottid in lookup: lookup[mrca] = lookup[ottid]
nwk = re.sub(r"([(,])(mrcaott\d+ott\d+|ott\d+)", lambda m: m.group(1) + lookup.get(m.group(2), m.group(2)), nwk)
nwk = re.sub(r"\)[^,;()]+", ")", nwk)
t = dendropy.Tree.get(data=nwk, schema="newick", preserve_underscores=True, rooting="force-rooted")
t.suppress_unifurcations()
tips = {l.taxon.label for l in t.leaf_node_iter()}
st.loc[st.sp.isin(cand.sp) & ~st.sp.isin(tips), "status"] = "lost_in_tree"

# graft missing species next to congeners
genus = dict(zip(st.sp, st.species.str.split().str[0]))
in_tree = {s for s in st.sp if s in tips}
for row in st[~st.sp.isin(in_tree) & (st.status != "duplicate_ott")].itertuples():
    cong = [s for s in in_tree if genus[s] == genus[row.sp]]
    if not cong:
        st.loc[st.sp == row.sp, "status"] += ";dropped"; continue
    node = t.mrca(taxon_labels=cong) if len(cong) > 1 else t.find_node_with_taxon_label(cong[0])
    if node is None:
        st.loc[st.sp == row.sp, "status"] += ";dropped"; continue
    tax = dendropy.Taxon(label=row.sp); t.taxon_namespace.add_taxon(tax)
    parent = node.parent_node
    new = dendropy.Node(); new.add_child(dendropy.Node(taxon=tax))
    if parent is None:
        t.seed_node = new; new.add_child(node)
    else:
        parent.remove_child(node); parent.add_child(new); new.add_child(node)
    st.loc[st.sp == row.sp, "status"] += ";grafted"
    in_tree.add(row.sp)
t.suppress_unifurcations()
(PROC / "species_tree_expanded.nwk").write_text(t.as_string(schema="newick").strip().replace("[&R] ", ""))
st[["sp", "species", "source", "ott_id", "status"]].to_csv(PROC / "species_tree_expanded_log.tsv", sep="\t", index=False)
print(st.status.replace("", "ok").value_counts().to_string())
print(len(t.leaf_nodes()), "tips")
