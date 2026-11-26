"""Download 1KP (OneKP capstone, 2019) orthogroup-matched peptides for every sample.

Source: CyVerse Data Commons, oneKP_capstone_2019/transcript_assemblies, via the
public Terrain API (anonymous; data.cyverse.org's web endpoint is bot-gated).
Per sample: <code>.faa.tar.bz2 = translated transcripts that matched a gene family
(31-genome OrthoFinder families, BLAST e < 1e-10), one FAA file per family.
Writes data/raw/onekp/<code>.faa.gz (all peptides, family id in header) and samples.tsv.
Usage: python3 08_fetch_onekp.py [n_workers]
"""
import io, sys, tarfile, gzip, pathlib, time, requests
from concurrent.futures import ThreadPoolExecutor

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "data/raw/onekp"; OUT.mkdir(parents=True, exist_ok=True)
BASE = "/iplant/home/shared/commons_repo/curated/oneKP_capstone_2019/transcript_assemblies"
API = "https://de.cyverse.org/terrain"
NW = int(sys.argv[1]) if len(sys.argv) > 1 else 6


def get(url, **kw):
    for i in range(5):
        try:
            r = requests.get(url, timeout=300, **kw); r.raise_for_status(); return r
        except Exception:
            if i == 4: raise
            time.sleep(5 * (i + 1))


def list_dir(path):
    out, off = [], 0
    while True:
        d = get(f"{API}/filesystem/paged-directory", params={"path": path, "limit": 500, "offset": off}).json()
        out += [(f["label"], "D") for f in d["folders"]] + [(f["label"], "F") for f in d["files"]]
        off += 500
        if off >= d["total"]: return out


def fetch(folder):
    code = folder.split("-")[0]
    dest = OUT / f"{code}.faa.gz"
    if dest.exists(): return folder, "cached"
    try:
        files = [n for n, t in list_dir(f"{BASE}/{folder}") if t == "F"]
        faa = next((f for f in files if f.endswith(".faa.tar.bz2")), None)
        if faa is None: return folder, "no faa"
        raw = get(f"{API}/fileio/download", params={"path": f"{BASE}/{folder}/{faa}"}).content
        tmp = dest.with_suffix(".tmp")
        with tarfile.open(fileobj=io.BytesIO(raw), mode="r:bz2") as tf, gzip.open(tmp, "wt") as out:
            for m in tf.getmembers():
                if not m.isfile(): continue
                og = m.name.rsplit("/", 1)[-1].split(".")[1]
                for line in tf.extractfile(m).read().decode().splitlines():
                    out.write(f"{line.split()[0]} og={og}\n" if line.startswith(">") else line + "\n")
        tmp.rename(dest)
        return folder, "ok"
    except Exception as e:
        return folder, f"error {e}"


if __name__ == "__main__":
    folders = [n for n, t in list_dir(BASE) if t == "D"]
    (OUT / "samples.tsv").write_text("folder\tcode\tlabel\n" + "".join(
        f"{f}\t{f.split('-')[0]}\t{f.split('-', 1)[1] if '-' in f else ''}\n" for f in folders))
    print(len(folders), "sample folders", flush=True)
    st = {}
    with ThreadPoolExecutor(NW) as ex:
        for i, (f, s) in enumerate(ex.map(fetch, folders)):
            st[s.split()[0]] = st.get(s.split()[0], 0) + 1
            if s.startswith("error") or s == "no faa": print(f, s[:200], flush=True)
            if i % 100 == 0: print(i, st, flush=True)
    print("finished", st)
