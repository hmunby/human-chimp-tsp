import json, time, urllib.request, urllib.parse, sys

BASE = "https://www.ebi.ac.uk/eqtl/api/v2"
DATASETS = [("minor_salivary_gland", "QTD000276", 144),
            ("skin_sun_exposed",     "QTD000316", 602)]

# The five TSPs of the MUC7 region (4_71210148_71210970 in hg19), GRCh38 ids.
# From ../../../4_SINGER_GenealogicalReconstruction_TSP/4_TSPRegions/
#      3_GeneAssignment/variant_gene_assignments.pooled.tsv
VARIANTS = [
    ("chr4_70344431_A_G", "rs34735123"),
    ("chr4_70344481_C_T", "rs13101613"),
    ("chr4_70344604_A_T", "rs67708234"),
    ("chr4_70344619_T_C", "rs13134277"),
    ("chr4_70345253_A_G", "rs13148087"),
]

def get(url):
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception as e:
            if attempt == 3:
                return {"__error__": str(e)}
            time.sleep(2 * (attempt + 1))

rows = []
for tissue, ds, n in DATASETS:
    for vid, rsid in VARIANTS:
        q = urllib.parse.urlencode({"variant": vid, "size": 1000})
        d = get(f"{BASE}/datasets/{ds}/associations?{q}")
        if isinstance(d, dict):
            rows.append({"tissue": tissue, "dataset": ds, "n": n, "variant": vid,
                         "rsid": rsid, "gene_id": None,
                         "note": d.get("message") or d.get("__error__") or str(d)[:200]})
        else:
            for a in d:
                rows.append({"tissue": tissue, "dataset": ds, "n": n, "variant": vid,
                             "rsid": rsid, "gene_id": a.get("gene_id"),
                             "pvalue": a.get("pvalue"), "beta": a.get("beta"),
                             "se": a.get("se"), "maf": a.get("maf"),
                             "ac": a.get("ac"), "an": a.get("an"),
                             "median_tpm": a.get("median_tpm"),
                             "pos": a.get("position"), "note": ""})
        time.sleep(0.3)

json.dump(rows, open("assoc_rows.json", "w"), indent=1)
print("rows:", len(rows))
